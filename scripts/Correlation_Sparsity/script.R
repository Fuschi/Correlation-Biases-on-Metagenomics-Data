library(tidyverse)
library(mvtnorm)
library(Matrix)
library(foreach)
library(doSNOW)
library(parallel)
library(ToyModel)


## Different Densities
#------------------------------------------------------------------------------#

# Set up the parallel backend with doSNOW
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

# Define the parameters
Dimension <- 500
Connected_Taxa_Sequence <- seq(200, 0, -1)
Correlation_Sequence <- seq(.5,.9,.1)

# Total possible edges
total_possible_edges <- Dimension * (Dimension - 1) / 2

# Set up progress bar
total_iterations <- length(Correlation_Sequence) * length(Connected_Taxa_Sequence)
pb <- txtProgressBar(max = total_iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Perform parallel computation with progress tracking
results_block <-  foreach(corr_val = Correlation_Sequence, .combine = "rbind") %:%
  foreach(n_taxa_connected = Connected_Taxa_Sequence, .combine = "rbind", .packages = c("mvtnorm", "Matrix", "ToyModel"), .options.snow = opts) %dopar% {
    
    # Start with a diagonal matrix of ones
    corM <- diag(1, Dimension, Dimension)
    
    # Get the current upper triangle values excluding the diagonal
    corM[1:n_taxa_connected, 1:n_taxa_connected] <- corr_val
    
    # Number of edges currently connected
    num_edges <- sum(corM[lower.tri(corM)] != 0)
    edge_density <- num_edges / total_possible_edges
    
    # Ensure the matrix is positive semidefinite
    corM_PD <- nearPD(corM, corr = TRUE, keepDiag = TRUE)$mat
    corM_PD <- as.matrix(corM_PD)
    
    toy <- toy_model(n=10^4, cor=corM_PD, M=1,
                     qdist=qnorm, 
                     param=c(mean=0, sd=1),
                     method="pearson",
                     force.positive=TRUE)
    
    data.frame("d"=Dimension, "taxa_connected" = n_taxa_connected,
               "edge_density" = edge_density,
               "correlation_value" = corr_val,
               "ERR_CLR"=mean(abs(toy$cor_NorTA - toy$cor_CLR)))
  }

# Close progress bar and stop the cluster
close(pb)
stopCluster(cl)

saveRDS(results_block, "results_block.rds")

## Gaussian Distributed Correlation Values
#------------------------------------------------------------------------------#

# Set up the parallel backend with doSNOW
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

# Set up progress bar
total_iterations <- 100
pb <- txtProgressBar(max = total_iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## Gaussian Random Values
n <- Dimension
results_gauss <- foreach(i = 1:100, .combine = "rbind", .packages = c("mvtnorm", "Matrix", "ToyModel"), .options.snow = opts) %dopar% {
  
  # Step 1: Generate a random matrix with Gaussian entries
  random_matrix <- matrix(rnorm(n*n), n, n)
  
  # Step 2: Make the matrix symmetric
  symmetric_matrix <- (random_matrix + t(random_matrix)) / 2
  
  # Step 3: Standardize the diagonal to make it a correlation matrix
  diag(symmetric_matrix) <- 1
  
  # Step 4: Ensure the matrix is positive semidefinite
  eigen_values <- eigen(symmetric_matrix)
  # Set any negative eigenvalues to a small positive value (e.g., 1e-10)
  eigen_values$values[eigen_values$values < 0] <- 1e-10
  # Reconstruct the matrix
  correlation_matrix <- eigen_values$vectors %*% diag(eigen_values$values) %*% t(eigen_values$vectors)
  
  # Verify the matrix is now a valid correlation matrix
  correlation_matrix <- nearPD(correlation_matrix, corr = TRUE, keepDiag = TRUE)$mat
  correlation_matrix <- as.matrix(correlation_matrix)
  
  toy <- toy_model(n=10^4, cor=correlation_matrix, M=1,
                   qdist=qnorm, 
                   param=c(mean=0, sd=1),
                   method="pearson",
                   force.positive=TRUE)
  
  data.frame("d" = n, "i" = i,
             "ERR_CLR" = mean(abs(toy$cor_NorTA - toy$cor_CLR)))
  
}

saveRDS(results_gauss, "results_gauss.rds")


# Create the plot with improved aesthetics
p_gauss <- results_gauss %>%
  ggplot(aes(x = factor(1), y = ERR_CLR)) +
  geom_violin(fill = "#69b3a2", color = "#1b4f4a", alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", color = "#2a2a2a", alpha = 0.9, outlier.shape = NA) +  # Hiding outliers as jitter will show them
  geom_jitter(width = 0.15, color = "#fc9272", alpha = 0.6, size = 2) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(y = "MAE") +
  geom_hline(yintercept = median(results_gauss$ERR_CLR), linetype = "dotdash", color = "darkgoldenrod1", linewidth = .75)
  #annotate("text", x = .2, y = Inf, label = "A", hjust = -0.1, vjust = 1.1, size = 6, color = "black")
p_gauss

p_block <- results_block %>%
  filter(edge_density <=.10) %>%
  filter(correlation_value %in% c(0.5, 0.7, 0.9)) %>%
  mutate(edge_density = 100*edge_density) %>%
  ggplot(aes(x = edge_density, y = ERR_CLR, color = as.factor(correlation_value))) + 
  geom_point() +
  ylab("MAE") +
  xlab("Correlation Density") +
  theme_bw() +
  theme(legend.position = "right") +
  labs(color = "Correlation Value") +
  xlab("Edge Density") +
  #annotate("text", x = 1, y = Inf, label = "B", hjust = 1.1, vjust = 1.1, size = 6, color = "black") +
  scale_x_continuous(labels = scales::label_percent(scale = 1)) +
  geom_hline(yintercept = median(results_gauss$ERR_CLR), linetype = "dotdash", color = "darkgoldenrod1", linewidth = 1.5) +
  annotate("text", x = 10, y = median(results_gauss$ERR_CLR), 
           label = sprintf("Median MAE from Gaussian-weighted\nCorrelation Matrices ≈ %.0e", median(results_gauss$ERR_CLR)), 
           hjust = 1.1, vjust = -.1, size = 4, color = "darkgoldenrod4") +
  labs(y = "MAE", x = "Edge Density (%)", color = "Correlation Value")
p_block



pall <- ggpubr::ggarrange(
  
  p_gauss, p_block, widths = c(.3, .7), labels = c("A", "B")
  
)

  
png(filename = "../Plots/correlation_density_biases.png", width = 2800, height = 1500, res = 300)  
pall
dev.off()
  
  
  
  
  
  
