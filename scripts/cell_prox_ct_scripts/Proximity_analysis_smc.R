## proximity code clean
library(sf)
library(dplyr)
library(data.table)
library(doParallel)
library(tidyr)
library(splancs)

filter <- dplyr::filter
select <- dplyr::select
pull <- dplyr::pull

## load functions 
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

## calculate distance function
calc_d <- function(pt1, pt2) {
  d <- sqrt( (( pt2[2] - pt1[2] )^2 ) + (( pt2[1] - pt1[1] )^2 ) )
  return(d)
}

## calculate angle
getDir <- function(origin, pt2){
  pt2 <- pt2 - origin
  deg <- abs(atan(pt2[2]/pt2[1]) * (180/pi))
  if(pt2[1] > 0 & pt2[2] > 0){deg <- deg
  } else if(pt2[1] < 0 & pt2[2] > 0){ deg <- 180-deg
  } else if(pt2[1] < 0 & pt2[2] < 0){ deg <- 180 + deg
  } else if(pt2[1] > 0 & pt2[2] < 0){ deg <- 360 - deg}
  return(deg)
}

# assign angle
assign.angle <- function(degrees, range){
  require(data.table)
  a=data.table(degrees=degrees)
  a[,merge:=degrees]
  
  b=data.table(range=range)
  b[,merge:=range]
  
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  Merge_a_b=b[a,roll='nearest']
  return(Merge_a_b)
}

anchorPoint <- function(degrees, anchor){
  object <- (degrees + (360 - anchor)) 
  object <- ifelse(object > 360, object - 360, object)
  return(object)
}


get_neighbors <- function(cell_x_coords, cell_y_coords, r, cell_coords_df, cell_id){
  
  all_cell_coord <- data.frame(x = cell_coords_df$x_centroid, y = cell_coords_df$y_centroid)
  circ_coord <- circleFun(c(cell_x_coords, cell_y_coords), r, 1000)
  circ_coord <- circ_coord[!duplicated(circ_coord),]
  # get cells within radius
  tmp <- circ_coord[,c(1:2)]
  io <- splancs::inout(all_cell_coord, tmp)
  
  if(mean(io) == 0){message("no proximal cells found"); return(NULL)}
  extract_cells <- as.data.frame(cell_coords_df)[which(io == TRUE),] ## get cell within radius
  if(nrow(extract_cells) == 1){message("no proximal cells found"); return(NULL)}
  
  # calculate distance and degrees
  ## Calculate distance between cells
  tmp <- extract_cells[,c('x_centroid', 'y_centroid')]
  tmp <- rbind(tmp[which(rownames(tmp) %in% cell_id),], tmp)
  tmp <- tmp[!duplicated(tmp),]
  
  toCheck <- combn(rownames(tmp), 2, simplify = FALSE)
  names(toCheck) <-
    sapply(toCheck, paste, collapse = " - ")
  toCheck <- toCheck[grep(pattern = cell_id, toCheck)]
  
  ## calculate distance
  tmp.d <- sapply(toCheck, function(j){
    calc_d(tmp[j[1],c(1,2)], tmp[j[2],c(1,2)]) })
  names(tmp.d) <- gsub('\\.y_centroid', '', names(tmp.d))
  tmp.d <- as.data.frame(do.call(rbind, tmp.d))
  tmp.d$cell_id <- rownames(tmp.d)
  tmp.d <- tidyr::separate(tmp.d, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
  rownames(tmp.d) <- NULL
  
  ## calculate degrees
  tmp.degree <- sapply(toCheck, function(j){
    getDir(c(cell_x, cell_y), tmp[j[2],c(1,2)]) })
  names(tmp.degree) <- gsub('\\.y_centroid', '', names(tmp.degree))
  tmp.degree <- as.data.frame(do.call(rbind, tmp.degree))
  tmp.degree$cell_id <- rownames(tmp.degree)
  tmp.degree <- tidyr::separate(tmp.degree, cell_id, sep = "\\ - ", c('cell.a', 'cell.b'))
  rownames(tmp.degree) <- NULL
  
  ## add degrees to distance table
  tmp.d$degree <- tmp.degree$V1[match(tmp.d$cell.b, tmp.degree$cell.b)]
  tmp.d$celltypeA <- extract_cells$CT_final[match(tmp.d$cell.a, rownames(extract_cells))]
  tmp.d$celltypeB <- extract_cells$CT_final[match(tmp.d$cell.b, rownames(extract_cells))]
  
  range <- seq(60, 360, 60) -30
  Merge_a_b <- assign.angle(degrees = tmp.d$degree, range = range)
  tmp.d$angle <- Merge_a_b$range[match(tmp.d$degree, Merge_a_b$degrees)]
  tmp.d <- tmp.d[!duplicated(tmp.d),]
  return(tmp.d)
}

cell_type_id <- 'SMC'
cell_type_id_out <- gsub("\\/", "\\_", cell_type_id)
print(cell_type_id)

## load data and define code parameters
inDir <- '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/input/'
outDir <- '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output/'

# read object metadata
x <- read.csv(file.path(inDir, "run2_final_metadata.csv"), row.names = 1)
x[1:5,1:5]
sample_ids <- x %>% 
  dplyr::pull(Sample) %>% as.character %>% unique() 

# set up enviorment to run in parallel
nCores <- 24
cl <- makeCluster(nCores)
registerDoParallel(cores=nCores)

# Export functions and packages to workers
cell_id = 0
cell_x = 0
cell_y = 0
clusterExport(cl = cl, 
              varlist = c("circleFun", "calc_d", "getDir", "assign.angle", "anchorPoint", 
                          "get_neighbors", "cell_id", "cell_x", "cell_y"),
              envir = environment())
clusterEvalQ(cl = cl, {
  library(sf)
  library(dplyr)
  library(data.table)
  library(tidyr)
})


## Loop over samples
proximal_cell_pop <- foreach (i=1:length(sample_ids)) %dopar% {
  
  sid <- sample_ids[i]; print(sid)
  obj.sid <- subset(x, Sample == sid)
  
  # 1. get cell coordinates
  meta <- obj.sid
  
  cells_coord <- meta %>% 
    filter(CT_final == cell_type_id)
  if(nrow(cells_coord) == 0){return(NULL)} ## skip is sample does not have the right cells
  
  cell_counts <- c()
  r <- 30 ## set radius ~30 is 3 cells 
  print(r)
  for(j in 1:nrow(cells_coord)){
    print(j)
    cell_id <- rownames(cells_coord)[j]
    split_cell.id <- strsplit(cell_id, ";") %>% unlist()
    cell_x <- cells_coord$x_centroid[j]
    cell_y <- cells_coord$y_centroid[j]
    
    tmp.d <- get_neighbors(cell_x_coords = cell_x,
                           cell_y_coords = cell_y,
                           r = r,
                           cell_coords_df = meta,
                           cell_id = cell_id)
    if(is.null(tmp.d)){next}
    tmp.d$cell.a <- cell_id
    tmp.d <- tmp.d[!duplicated(tmp.d),]
    
    # index rows by distance + angle
    tmp.d <- tmp.d %>%
      group_by(angle) %>%
      mutate(idx = rank(V1, ties.method = "min")) %>%
      ungroup()
    
    cell_counts <- rbind(cell_counts, tmp.d)
  }
  cell_counts$sid <- sid
  return(cell_counts)
}

proximal_cell_pop_clean <- proximal_cell_pop[
  sapply(proximal_cell_pop, function(x) is.data.frame(x) && !is.null(x))
]

proximal_cell_pop_clean <- do.call("rbind", proximal_cell_pop_clean)
stopImplicitCluster() 
out_path <- '/scratch/smallapragada/run2/proximity_results_allsamples_allCTs_2025_04_18/output'
out_fname <- paste0(cell_type_id_out, '_cell_dist_adjdegree.rds')
saveRDS(proximal_cell_pop_clean, file.path(out_path, out_fname)) 
