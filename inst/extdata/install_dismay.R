



devtools::install_github("tpq/propr") # propr: https://www.nature.com/articles/s41598-017-16520-0
BiocManager::install("preprocessCore")
BiocManager::install("impute")
BiocManager::install("GO.db")

devtools::install_github("skinnider/dismay")

dismay::kendall_zi()
