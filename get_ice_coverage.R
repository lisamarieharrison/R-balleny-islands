#gets ice coverage for Balleny Islands data
#code from https://github.com/AustralianAntarcticDivision/raadtools


devtools::install_github("AustralianAntarcticDivision/raadtools")
devtools::install_github("AustralianAntarcticDataCentre/raadsync")

cfg <- read_repo_config(system.file("extdata", "raad_repo_config.json", package= "raadsync"))
ice_index <- 1
cfg$do_sync <- seq(nrow(cfg)) == ice_index
## limit our data to only the last year
cfg$method_flags[1] <- paste0(cfg$method_flags[1], " --accept=\"*nt_201502*\"")
## specify local repository location
## NOTE: this should point to your permanent data collection - at whatever path you need
my_datadir <- normalizePath("~/Lisa/phd/Balleny Islands/remote data/sea ice", "/")
options(default.datadir = my_datadir)
cfg$local_file_root <- file.path(my_datadir, "data")

## 2. trigger raadsync to build the repository (this is set up in a cron-task ultimately)
er <- sync_repo(cfg, create_root = TRUE)
## build the file cache (many, many files for some data, so needs a system shortcut )
my_datadir <- getOption("default.datadir")
my_admindir <- file.path(my_datadir, 'admin', 'filelist')
if (!file.exists(my_admindir)) dir.create(my_admindir, recursive = TRUE)
fs <- list.files(file.path(my_datadir, 'data'), all = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)
fs <- gsub(pattern=".*sea ice/", replacement = "", x=fs)
save(fs, file = file.path(my_datadir, 'admin', 'filelist', 'allfiles2.Rdata'))
writeLines(fs, file.path(my_datadir, 'admin', 'filelist', 'allfiles2.txt'))
