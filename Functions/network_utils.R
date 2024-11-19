
local_file <- withr::local_tempfile()

#' @export
string_db_local <- function(file_path = local_file) {
  
  file_path |> 
    readr::read_delim(delim = " ") |> 
    # remove combined score, as we are calculating ourselves
    dplyr::select(-"combined_score")
}

#' @export
string_db_homo_sapiens <- function(
    version = "12.0", 
    species = 9606,
    links = "links.full",
    url.domain = "stringdb-downloads.org"
) {
  url <- glue::glue(
    .sep = "/",
    "https://{url.domain}",
    "download/protein.{links}.v{version}",
    "{species}.protein.{links}.v{version}.txt.gz"
  )
  if (!file.exists(local_file)) {
    curl::curl_download(url, local_file)
  }
  invisible()
}

#' @export
my_plot <- \(x) {
  lambda_min_count <- sum(glmnet::coef.glmnet(x, s = "lambda.min") != 0)
  lambda_1se_count <- sum(glmnet::coef.glmnet(x, s = "lambda.1se") != 0)
  title_str <- glue::glue("{attr(x, 'label')} (selected vars min: {lambda_min_count} ; 1se: {lambda_1se_count})")
  graphics::plot(x)
  graphics::title(sub = title_str)
}