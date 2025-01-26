#' A function to add accesssory data to a data-list for use in multiplex analysis
#'
#' These are internal functions.
#'
#' @param data A data list.
#' @return A data list with extra information about the structure of the dyadic or generalized reciprocity matrix.
#' @export
#'

build_multiplex_bindings_dr_multiplex = function(data){
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = seq(1,choose(data$N_responses,2))
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0

  # Build B matrix
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[lower.tri(bindings_mat_B)] = seq(1,choose(data$N_responses,2)) + choose(data$N_responses,2)
  bindings_mat_B[upper.tri(bindings_mat_B)] = 0
  bindings_mat_B = bindings_mat_B + t(bindings_mat_B)
  
  # Assemble full matrix
  bindings_mat = diag(NA, 2*data$N_responses)
  bindings_mat[upper.tri(bindings_mat)] = 0

  bindings_mat[1:data$N_responses, 1:data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, (1:data$N_responses)+data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, 1:data$N_responses] = bindings_mat_B
  
  # Compute lookup table
  dr_indices = which(bindings_mat > 0, arr.ind=TRUE)
  dr_indices = dr_indices[order(dr_indices[,1]),]
  dr_id = bindings_mat[dr_indices]
  
  # Compute other info 
  data$N_dr_bindings = 2*choose(data$N_responses,2)
  data$N_dr_params = data$N_dr_bindings + data$N_responses
  data$N_dr_indices = nrow(dr_indices)
  data$dr_indices = dr_indices
  data$dr_id = bindings_mat[dr_indices]

  print(bindings_mat)

  return(data)
 }


build_multiplex_bindings_dr_longitudinal = function(data){
  # Compute banding structure
  fill_set = c()
   for(k in 1:(data$N_responses-1)){
    fill_set = c(fill_set, seq(k, 1))
   }
  
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = fill_set
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0
  
  # Build B matrix
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[upper.tri(bindings_mat_B)] = fill_set + data$N_responses - 1
  bindings_mat_B[lower.tri(bindings_mat_B)] = 0
  bindings_mat_B = bindings_mat_B + t(bindings_mat_B)
  diag(bindings_mat_B) = 2*data$N_responses - 1

  # Assemble full matrix
  bindings_mat = diag(NA, 2*data$N_responses)
  bindings_mat[upper.tri(bindings_mat)] = 0

  bindings_mat[1:data$N_responses, 1:data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, (1:data$N_responses)+data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, 1:data$N_responses] = bindings_mat_B

  # Compute lookup table
  dr_indices = which(bindings_mat > 0, arr.ind=TRUE)
  dr_indices = dr_indices[order(dr_indices[,1]),]
  dr_id = bindings_mat[dr_indices]

  # Compute other info 
  data$N_dr_bindings = choose(2*data$N_responses,2) # Not sure about this one
  data$N_dr_params = 2*data$N_responses - 1
  data$N_dr_indices = nrow(dr_indices)
  data$dr_indices = dr_indices
  data$dr_id = bindings_mat[dr_indices]

  print(bindings_mat)

  return(data)
 }
 

 build_multiplex_bindings_sr_longitudinal = function(data){
  # Compute banding structure
  fill_set = c()
   for(k in 1:(data$N_responses-1)){
    fill_set = c(fill_set, seq(k, 1))
   }
  
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = fill_set
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0

  # Build C matrix
  bindings_mat_C = diag(NA, data$N_responses)
  bindings_mat_C[upper.tri(bindings_mat_C)] = fill_set + (data$N_responses - 1)
  bindings_mat_C[lower.tri(bindings_mat_C)] = 0
  
  # Build B matrice
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[upper.tri(bindings_mat_B)] = fill_set + 2*(data$N_responses - 1)
  bindings_mat_B[lower.tri(bindings_mat_B)] = 0

  bindings_mat_Bt = diag(0, data$N_responses)
  bindings_mat_Bt[upper.tri(bindings_mat_Bt)] = fill_set + 3*(data$N_responses - 1)
  bindings_mat_Bt[lower.tri(bindings_mat_Bt)] = 0

  bindings_mat_B = bindings_mat_B + t(bindings_mat_Bt)
  diag(bindings_mat_B) = 4*(data$N_responses - 1) + 1

  # Assemble full matrix
  bindings_mat = diag(NA, 2*data$N_responses)
  bindings_mat[upper.tri(bindings_mat)] = 0

  bindings_mat[1:data$N_responses, 1:data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, (1:data$N_responses)+data$N_responses] = t(bindings_mat_C)
  bindings_mat[(1:data$N_responses)+data$N_responses, 1:data$N_responses] = bindings_mat_B

  # Compute lookup table
  sr_indices = which(bindings_mat > 0, arr.ind=TRUE)
  sr_indices = sr_indices[order(sr_indices[,1]),]
  sr_id = bindings_mat[sr_indices]

  # Compute other info 
  data$N_sr_bindings = choose(2*data$N_responses,2) # Not sure about this one
  data$N_sr_params = 4*(data$N_responses - 1) + 1
  data$N_sr_indices = nrow(sr_indices)
  data$sr_indices = sr_indices
  data$sr_id = bindings_mat[sr_indices]

  print(bindings_mat)

  return(data)
 }


 build_multiplex_bindings_free = function(data){
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = 0
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0

  # Build B matrix
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[lower.tri(bindings_mat_B)] = 0
  bindings_mat_B[upper.tri(bindings_mat_B)] = 0
  bindings_mat_B = bindings_mat_B + t(bindings_mat_B)
  
  # Assemble full matrix
  bindings_mat = diag(NA, 2*data$N_responses)
  bindings_mat[upper.tri(bindings_mat)] = 0

  bindings_mat[1:data$N_responses, 1:data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, (1:data$N_responses)+data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, 1:data$N_responses] = bindings_mat_B
  
  # Compute lookup table
  dr_indices = which(bindings_mat > 0, arr.ind=TRUE)
  dr_indices = dr_indices[order(dr_indices[,1]),]
  dr_id = bindings_mat[dr_indices]
  
  # Compute other info 
  data$N_dr_bindings = 0
  data$N_dr_params = choose(2*data$N_responses,2)
  data$N_dr_indices = 0
  data$dr_indices = dr_indices
  data$dr_id = bindings_mat[dr_indices]

  print(bindings_mat)

  return(data)
 }
