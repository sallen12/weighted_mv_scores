################################################################################
# helper functions

# output checks for chaining function in threshold-weighted scores
check.chaining <- function(x1, x2) {
  
  if (is.vector(x1) && is.vector(x2)) {
    if (length(x1) != length(x2)) {
      stop("The chaining function changes the length of the input")
    }
  }else if (is.matrix(x1) && is.matrix(x2)) {
    if (any(dim(x1) != dim(x2))) {
      stop("The chaining function changes the dimensions of the input")
    }
  }else {
    warning("The chaining function changes the structure of the data")
  }
  
}

# output checks for weight function in vertically re-scaled scores
check.mv.weight <- function(dat, w) {
  if(is.vector(dat)){
    m <- 1
  }else{
    m <- ncol(dat)
  }
  if (is.null(w)) {
    w_final <- rep(1, m)
  } else {
    if (any(w < 0) || (length(w) != m) || !is.vector(w)){
      stop("w is of wrong format")
    }
    w_final <- w
  }
  return(w_final)
}

# output checks for weight function in vertically re-scaled scores
check.uv.weight <- function(dat, w) {
  if (is.matrix(dat)) {
    m <- dim(dat)
  } else if (identical(length(dat), 1L)) {
    m <- 1
  } else {
    m <- length(dat)
  }
  
  if (is.null(w)) {
    w_final <- array(1, m)
  } else {
    if (any(w < 0) || !is.vector(w) & is.vector(dat) || is.vector(w) & !is.vector(dat) || is.vector(dat) & (length(w) != m)){
      stop("w is of wrong format")
    }
    w_final <- w
  }

  return(w_final)
}

# output checks for weight function in variogram score
w_vs.helper <- function(dat, w_vs) {
  m <- ncol(dat)
  if (!is.null(w_vs)) {
    if (!is.matrix(w_vs)) {
      stop("'w_vs' is not a matrix ")
    }
    if (any(dim(w_vs) != d)) {
      stop("Dimensions of 'w_vs' do not fit")
    }
    if (any(w_vs < 0)) {
      stop("Weighting matrix 'w_vs' contains negative values")
    }
    if (!isSymmetric(w_vs)) {
      stop("Weighting matrix 'w_vs' is not symmetric")
    }
  }else{
    w_vs <- matrix(1, m, m)
  }
  return(w_vs)
}


################################################################################
# helper functions from scoringRules

# input checks for multivariate scoring rules
check.multivsample <- function(input) {
  input_isnumeric <- sapply(input, is.numeric)
  if (!all(input_isnumeric)) {
    stop(paste("Non-numeric input:", 
               paste(names(input)[!input_isnumeric], collapse=", ")))
  }
  if (!is.vector(input$y)) {
    stop("'y' is not a vector")
  } 
  if (!is.matrix(input$dat)) {
    stop("'dat' is not a matrix ")
  }
  if (length(input$y) != dim(input$dat)[1]) {
    stop("Dimensions of 'y' and 'dat' do not fit")
  }
}

# function to check or generate weights
w.helper.multiv <- function(dat, w){
  m <- ncol(dat)
  if (is.null(w)){
    w_final <- rep(1, m)/m
  } else {
    if (any(w < 0) || (length(w) != m) || !is.vector(w)){
      stop("w is of wrong format")
    } else {
      sum_w <- sum(w)
      if (abs(sum_w - 1) > 1e-5){
        message("Weights in w don't add to one - they have been re-scaled to one")
      }
      w_final <- w/sum_w
    }
  }
  return(w_final)
}

# function to check or generate centre
x0.helper.multiv <- function(dat, x0){
  d <- nrow(dat)
  if (is.null(x0)){
    x0_final <- rep(0, d)
  } else if (length(x0) == d){
    x0_final <- x0
  } else if (length(x0) == 1) {
    x0_final <- rep(x0, d)
  } else {
      stop("x0 is of wrong format")
    }
  return(x0_final)
}

# input checks for univariate sample functions
check_sample <- function(input) {
  checkNumeric(input)
  
  input_isvector <- sapply(input, is.vector)
  if (!all(input_isvector)) {
    stop(paste("Non-vector input:",
               paste(names(input)[!input_isvector], collapse=", ")))
  }
  
  input_lengths <- sapply(input, length)
  max_length <- max(input_lengths)
  ref_lengths <- c(y = 1L, dat = max_length, w = max_length, bw = 1L)
  if (!identical(input_lengths, ref_lengths[names(input)])) {
    ref_lengths2 <- c(y = "1", dat = "n", w = "n", bw = "1")
    stop(
      paste(
        "Incompatible input vector lengths.",
        sprintf("Lengths of (%s) should be (%s).",
                paste(names(input), collapse = ", "),
                paste(ref_lengths2[names(input)], collapse = ", ")),
        sprintf("Given lengths: %s", paste(input_lengths, collapse = ", ")),
        sep = "\n")
    )
  }
  
  if (!is.null(input$w)) {
    if (any(input$w < 0)) {
      stop("Weight parameter 'w' contains negative values.")
    }
  }
  if (!is.null(input$bw)) {
    if (input$bw < 0) {
      stop("Bandwidth parameter 'bw' is negative.")
    }
  }
}

check_sample2 <- function(input) {
  checkNumeric(input)
  
  input_isvector <- sapply(input[names(input) %in% c("y", "bw")], is.vector)
  input_ismatrix <- sapply(input[names(input) %in% c("dat", "w")], is.matrix)
  input_dim1 <- sapply(input, function(x) {
    if (is.vector(x)) length(x) else dim(x)[1L]
  })
  input_dim2 <- sapply(input[names(input) %in% c("dat", "w")], function(x) {
    if (is.matrix(x)) dim(x)[2L]
  })
  
  if (!all(input_isvector) ||
      !all(input_ismatrix) ||
      !identical(length(unique(input_dim1)), 1L) ||
      !identical(length(unique(input_dim2)), 1L)) {
    
    reference_formats <- c(
      y = "vector[1:n]",
      dat = "matrix[1:n, 1:m]",
      w = "matrix[1:n, 1:m]",
      bw = "vector[1:n]"
    )
    input_formats <- sapply(input, function(x) {
      if (is.vector(x)) {
        sprintf("vector[1:%i]", length(x))
      } else if (is.matrix(x)) {
        sprintf("matrix[1:%i, 1:%i]", dim(x)[1L], dim(x)[2L])
      } else {
        "unidentified"
      }
    })
    stop(
      paste(
        "Incompatible input objects.",
        sprintf("Expected input for (%s): %s.",
                paste(names(input), collapse = ", "),
                paste(reference_formats[names(input)], collapse = ", ")),
        sprintf("   Given input for (%s): %s.",
                paste(names(input), collapse = ", "),
                paste(input_formats, collapse = ", ")),
        sep = "\n"
      )
    )
  }
  
  if (!is.null(input$w)) {
    if (any(input$w < 0)) {
      stop("Weight parameter 'w' contains negative values.")
    }
  }
  if (!is.null(input$bw)) {
    if (input$bw < 0) {
      stop("Bandwidth parameter 'bw' is negative.")
    }
  }
}

checkNumeric <- function(input, infinite_exception = NULL) {
  input_numeric <- sapply(input, is.numeric)
  if (any(!input_numeric)) {
    stop(paste("Non-numeric input:",
               paste(names(input)[!input_numeric], collapse = ", ")))
  }
  input_NA <- sapply(input, anyNA)
  if (any(input_NA)) {
    stop(paste("Input with missing values:",
               paste(names(input)[input_NA], collapse = ", ")))
  }
  input <- input[!names(input) %in% infinite_exception]
  input_infinite <- sapply(input, function(x) any(is.infinite(x)))
  if (any(input_infinite)) {
    stop(paste("Input with infinite values:",
               paste(names(input)[input_infinite], collapse = ", ")))
  }
}
