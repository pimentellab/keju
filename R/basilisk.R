library(basilisk)

processing_env <- BasiliskEnvironment(envname="processing_keju",
    pkgname="keju",
    packages=c(
               "python=3.13.3",
               "numpy==2.2.5",
               "pandas==2.2.3", 
               "formulaic==1.1.1")
)

# .onLoad <- function() {
#     basilisk::configureBasiliskEnv()
# }