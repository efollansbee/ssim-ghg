find.indir <- function(settings.yaml="../../site_settings.yml") {
    library("yaml")
    settings <- read_yaml(settings.yaml)
    return(settings$global_paths$input_folder)
}

