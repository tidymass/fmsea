#' @title Show the base information of fmsea pacakge
#' @description Show the base information of fmsea pacakge.
#' @author Xiaotao Shen
#' \email{xiaotao.shen@@outlook.com}
#' @return A ASCII log of fmsea
#' @importFrom magrittr %>%
#' @importFrom crayon red yellow green bgRed
#' @importFrom stringr str_detect str_extract str_extract_all
#' @importFrom stringr str_replace_all str_replace str_trim str_c str_count
#' @importFrom readr cols read_csv
#' @importFrom pbapply pblapply pboptions
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr mutate filter everything select bind_rows left_join pull
#' @importFrom plyr dlply .
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam
#' @importFrom readxl read_xlsx read_xls
#' @importFrom purrr map map2
#' @importFrom ggplot2 aes ggplot geom_point geom_line geom_smooth theme annotate
#' @importFrom ggplot2 geom_abline theme_bw ggsave geom_segment xlim ylim labs
#' @importFrom ggplot2 element_line element_text
#' @importFrom MSnbase readMSData
#' @importFrom ProtGenerics spectra
#' @importFrom masstools ms2_match get_os mz_rt_match read_mzxml read_mgf
#' @importFrom masstools get_spectra_match_score
#' @importFrom stats lm loess predict
#' @importFrom plotly ggplotly
#' @importFrom future plan multisession
#' @importFrom furrr future_map2 future_map
#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom cli rule symbol cat_line
#' @importFrom data.table rbindlist
#' @importFrom progress progress_bar
#' @importFrom tidyr pivot_longer
#' @import RColorBrewer
#' @import utils
#' @import ggplot2
#' @import methods
#' @import graphics
#' @import grDevices
#' @import utils
#' @importClassesFrom massdataset mass_dataset
#' @export
#' @examples
#' fmsea_logo()

fmsea_logo <- function() {
  cat(crayon::green("Thank you for using fmsea!\n"))
  message(crayon::green("Version", fmsea_version, "(", update_date, ')\n'))
  cat(crayon::green("More information: google tidymass fmsea.\n"))
  cat(crayon::green(
    c(
      "   __ __  __  _____ ______          ",
      "  / _|  \\/  |/ ____|  ____|   /\\    ",
      " | |_| \\  / | (___ | |__     /  \\   ",
      " |  _| |\\/| |\\___ \\|  __|   / /\\ \\  ",
      " | | | |  | |____) | |____ / ____ \\ ",
      " |_| |_|  |_|_____/|______/_/    \\_\\",
      "                                    ",
      "                                    "
    )
    
  ), sep = "\n")
}

fmsea_version <- utils::packageVersion(pkg = "fmsea")
update_date <- as.character(Sys.time())

# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# fmsea_logo <-
#   c("                _    _____  ___ ", " _ __ ___   ___| |_  \\_   \\/   \\",
#     "| '_ ` _ \\ / _ \\ __|  / /\\/ /\\ /", "| | | | | |  __/ |_/\\/ /_/ /_// ",
#     "|_| |_| |_|\\___|\\__\\____/___,'  ", "                                "
#   )
# cat(fmsea_logo, sep = "\n")
