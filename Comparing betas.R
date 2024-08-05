library(dplyr)

df1K <- read.csv("../../nfs_share/Milner/Psupertime/With denoising/1K/images/1KpsupertimeResults.csv")
df2K <- read.csv("../../nfs_share/Milner/Psupertime/With denoising/2K/images/2KpsupertimeResults.csv")

min_beta <- 0.1
df1K <- df1K %>%
  filter(abs_beta >= min_beta)

min_beta <- 0.1
df2K <- df2K %>%
  filter(abs_beta >= min_beta)


common_genes <- intersect(df1K$symbol, df2K$symbol)

unique_1K_genes <- setdiff(df1K$symbol, df2K$symbol)

unique_2K_genes <- setdiff(df2K$symbol, df1K$symbol)


# Genes in both lists with their betas
common_df <- df1K %>%
  filter(symbol %in% common_genes) %>%
  select(symbol, beta_1K = beta) %>%
  inner_join(df2K %>% filter(symbol %in% common_genes) %>% select(symbol, beta_2K = beta), by = "symbol")

# Genes only in 1K
unique_1K_df <- df1K %>%
  filter(symbol %in% unique_1K_genes) %>%
  select(symbol, beta_1K = beta)

# Genes only in 2K
unique_2K_df <- df2K %>%
  filter(symbol %in% unique_2K_genes) %>%
  select(symbol, beta_2K = beta)


write.csv(common_df, file = "../../nfs_share/Milner/Psupertime/With denoising/Comparison/CommonGeneswithBetas.csv")
write.csv(unique_1K_df, file = "../../nfs_share/Milner/Psupertime/With denoising/Comparison/GenesUniqueto1KWithBetas.csv")
write.csv(unique_2K_df, file = "../../nfs_share/Milner/Psupertime/With denoising/Comparison/GenesUniqueto2KWithBetas.csv")





