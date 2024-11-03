#######################
####mavd_source####

# Function to read a PCL file
read.pcl <- function(filename, na.type = "", Nrows = -1, Comment.char = "", ...) {
  # Read the table from the PCL file, setting the first column as row names
  x.df <- read.table(paste(filename, "pcl", sep = "."), header = TRUE, sep = "\t",
                     quote = "\"", na.strings = na.type, skip = 0, nrows = Nrows, comment.char = Comment.char, ...);
  rownames(x.df) <- x.df[[1]];  # Set row names
  return(x.df);
};

# Function to write a data frame to a PCL file
write.pcl <- function(df, dataname, fileaddress = "") {
  # Construct the file path for the output PCL file
  dir.address <- paste(fileaddress, dataname, ".pcl", sep = "");
  # Write the data frame to a PCL file
  X <- write.table(df, file = dir.address,
                   append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "", row.names = FALSE, col.names = TRUE);
  return(X);
};

##########################
#### Various Functions ####
# Function to select elements from a vector based on indices
select.v <- function(x, indx) { 
  y <- x[indx]; 
  return(y); 
};

# Function to match and filter rows of two data frames by their row names
meshRows.mavd <- function(df1, df2) {
  rn <- rownames(df1)[{rownames(df1) %in% rownames(df2)}];  # Find common row names
  ndf1 <- df1[rn,]; 
  ndf2 <- df2[rn,]; 
  return(list(ndf1, ndf2));  # Return the filtered data frames
};

# Function to calculate the Euclidean norm of a vector
normvec <- function(vec) { 
  norm <- sqrt(sum(vec * vec)); 
  return(norm); 
};

# Function to convert a matrix to PCL format
mat2pcl <- function(mat, tag) { 
  mat.df <- as.data.frame(mat); 
  mat.pcl <- cbind(tag, mat.df);  # Combine tag with the data frame
  return(mat.pcl); 
};

#### PCA and Data Normalization Functions ####
# Function to create a diagonal matrix or convert a vector to a matrix
true.diag <- function(vec) { 
  ifelse(length(vec) > 1.5, y <- diag(vec), y <- as.matrix(vec)); 
  return(y); 
};

# Function to modify a vector by replacing values after a specified index
plugnew.vec <- function(vec, dimchange, newval = 0) { 
  vecnew <- vec; 
  vecnew[-(1:dimchange)] <- rep(newval, length(vec) - dimchange); 
  return(vecnew); 
};

# Function to perform PCA and return a list of matrices for each principal component
PCA.collapse.mat <- function(mat) {
  mat.svd <- svd(mat);  # Perform SVD on the input matrix
  U.mat <- mat.svd$u; 
  V.mat <- mat.svd$v; 
  D.mat <- list();  # List to hold diagonal matrices
  
  # Generate diagonal matrices
  for(j in 1:length(mat.svd$d)) {
    D.mat[[j]] <- true.diag(plugnew.vec(vec = mat.svd$d, dimchange = j));
  }
  
  N <- list();  # List to hold reduced dimension matrices
  
  # Calculate reduced dimension matrices
  for(j in 1:ncol(U.mat)) { 
    N[[j]] <- U.mat %*% D.mat[[j]] %*% t(V.mat); 
  }
  # Preserve original matrix attributes
  NN <- lapply(N, change.attributes, new.atr = attributes(mat)); 
  return(NN); 
};

# Function to collapse a matrix using specified number of principal components
PCA.1collapse.mat <- function(mat.svd, indx) { 
  U.mat <- mat.svd$u; 
  V.mat <- mat.svd$v;
  D.mat <- true.diag(plugnew.vec(vec = mat.svd$d, dimchange = indx));  # Create a diagonal matrix for specified index
  N <- U.mat %*% D.mat %*% t(V.mat);  # Generate the reduced matrix
  return(N); 
};

# Function to extract the first principal component as a vector
PC1.col.vec <- function(mat) { 
  svd1 <- svd(mat, nu = 1, nv = 1);  # Perform SVD, extracting one component
  u.vec <- as.vector(svd1$u, mode = "numeric");
  lambda <- svd1$d[[1]]; 
  v1 <- svd1$v[1,]; 
  N <- lambda * v1 * u.vec; 
  names(N) <- rownames(mat); 
  return(N); 
};

# Function to match rows of two data frames and normalize their values
meshRows.norm.mavd <- function(df1, df2, meanFUN = mean, ...) {
  # Call meshRows.mavd to match rows
  junk <- meshRows.mavd(df1, df2); 
  Df1 <- junk[[1]]; 
  Df2 <- junk[[2]]; 
  rm(junk);  # Remove intermediate variable
  
  # Extract matrices and calculate mean norms
  mat1 <- as.matrix(Df1); 
  mat2 <- as.matrix(Df2); 
  junk1 <- mean(apply(mat1, 2, normvec)); 
  junk2 <- mean(apply(mat2, 2, normvec)); 
  junk <- meanFUN(c(junk1, junk2), ...);  # Calculate overall mean norm
  
  # Normalize the data frames
  Df1 <- apply(Df1, 2, function(x) x / junk);
  Df2 <- apply(Df2, 2, function(x) x / junk);
  return(list(Df1, Df2)); 
};

##############################
##### Linear Model Functions ####
# Function to flatten the input matrix using linear regression
flat.mavd <- function(mat) {
  f1 <- function(x, m) { 
    nx <- fitted(lm(x ~ (m - 1))); 
    return(nx); 
  };
  f2 <- function(j, m) { 
    v <- m[, j]; 
    mm <- m[, -(j)]; 
    nv <- f1(v, mm); 
    return(nv); 
  };
  bigindx <- rbind(1:ncol(mat)); 
  colnames(bigindx) <- colnames(mat); 
  matt <- apply(bigindx, 2, f2, m = mat); 
  return(matt); 
};

# Function to calculate Wold's T² statistic for evaluating component importance
wold.mavd <- function(mat) {
  f1 <- function(svd.list, l) {
    n <- nrow(svd.list$u); 
    k <- ncol(svd.list$u); 
    lam.square <- {svd.list$d}^2; 
    junk <- {{lam.square[[l]]} / {sum(lam.square[-(1:l)])}} * {{n - l - 1} * {k - l} / {n + k - (2 * l)}};
    return(junk); 
  };
  svd.list <- svd(mat); 
  dim.list <- as.list((1:(length(svd.list$d) - 1))); 
  junk <- lapply(dim.list, f1, svd.list = svd.list); 
  junk <- unlist(junk); 
  return(junk); 
};

# Function to plot Wold's T² statistics
plot.wold.mavd <- function(x, x.Lbound = 1, x.Ubound = length(x), main.extra = "") { 
  plot.range <- (x.Lbound:x.Ubound); 
  z <- x[plot.range];
  # Create a line plot for Wold's T² statistic
  plot(plot.range, z, type = "l", lty = "solid", xlab = "Dimension of PC space = l", ylab = "W(l)",
       main = paste("Wold invariant", main.extra), col = "blue", log = "y"); 
  return(y); 
};

# Function to perform PCA and retain specified principal components
pca.mavd <- function(mat, j) {
  td <- function(vec) { 
    if(length(vec) > 1.5) 
      y = diag(vec) 
    else 
      y = as.matrix(vec); 
    return(y); 
  };
  pn <- function(vec) { 
    vn <- vec; 
    vn[-(1:j)] <- rep(0, length(vec) - j); 
    return(vn); 
  };
  mat.svd <- svd(mat); 
  U.mat <- mat.svd$u; 
  V.mat <- mat.svd$v; 
  D.mat <- td(pn(vec = mat.svd$d)); 
  matt <- U.mat %*% D.mat %*% t(V.mat); 
  attributes(matt) <- attributes(mat);  # Preserve original attributes
  return(matt); 
};

# Function for leave-one-out PCA
leave1out.mavd <- function(mat) {
  f1 <- function(vec) { 
    return(PC1.col.vec(mat[, vec, drop = FALSE])); 
  };
  bigindx <- rbind(1:ncol(mat)); 
  colnames(bigindx) <- colnames(mat); 
  matt <- apply(bigindx, 2, f1);  # Apply function to each column
  return(matt); 
};

# Function to normalize columns of a matrix against a linear model
normal.mavd <- function(mat, model) {
  coeff <- coefficients(lm(mat ~ model)); 
  normalized <- as.matrix(model) %*% coeff; 
  return(mat - normalized); 
};

# Function to extract coefficients from a linear regression model
normal.coefficients.mavd <- function(mat, model) {
  return(coefficients(lm(mat ~ model))); 
};

# Function to calculate residuals from linear regression
disease.mavd <- function(mat1, mat2) {
  N <- normal.mavd(mat1, mat2);  # Normalize matrix 1 against matrix 2
  return(N); 
};

##############################
### Main Analysis Functions ###
# Function to process normal and tumor data for analysis
mavd_part1 <- function(df.normal, df.tumor) {
  # Normalize and mesh rows of the input data
  junk <- meshRows.norm.mavd(df.normal, df.tumor); 
  df1 <- junk[[1]]; 
  df2 <- junk[[2]]; 
  
  # Perform PCA on the normalized data
  mat1 <- PCA.collapse.mat(as.matrix(df1)); 
  mat2 <- PCA.collapse.mat(as.matrix(df2)); 
  wold_results <- wold.mavd(as.matrix(df1));  # Calculate Wold's T²
  return(list(mat1 = mat1, mat2 = mat2, wold = wold_results)); 
};

# Function to continue analysis from part 1
mavd_part2 <- function(df.normal, df.tumor) {
  junk <- mavd_part1(df.normal, df.tumor); 
  mat1 <- junk$mat1; 
  mat2 <- junk$mat2; 
  # Perform further PCA and leave-one-out analysis
  leave1out_results <- leave1out.mavd(mat1); 
  return(list(leave1out = leave1out_results)); 
}
