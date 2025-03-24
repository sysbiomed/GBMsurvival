
local_file <- withr::local_tempfile()

#' @export
memoise_cache <- cachem::cache_disk(here::here("run-cache"))

#' @export
string_db_local_cache <- memoise::memoise(
  function(file_path = local_file) {
    
    file_path |> 
      readr::read_delim(delim = " ") |> 
      # remove combined score, as we are calculating ourselves
      dplyr::select(-"combined_score")
  },
  cache = memoise_cache
)

#' @export
calculate.combined.score_cache <- memoise::memoise(
  glmSparseNet:::calculate.combined.score,
  cache = memoise_cache
)

#' @export
buildStringNetwork_cache <- memoise::memoise(
  glmSparseNet::buildStringNetwork,
  cache = memoise_cache
)

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