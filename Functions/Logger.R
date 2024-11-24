log_chunk_end_time <- function(..., log_file = "ignored", namespace = "chunk") {
  logger::log_level(CHUNK, ..., namespace = namespace)
}

# Execute on environment of `log_chunk_end_time`
local({
  pryr::mem_used

  purple <- crayon::combine_styles(crayon::bold, crayon::make_style("purple"))
  gray <- crayon::make_style("gray70")

  CHUNK <- structure(5L, level = "CHUNK", class = c("loglevel", "integer"))

  logger::log_threshold(CHUNK, namespace = "chunk")

  local_layout <- logger::layout_glue_generator(
    format = paste(
      "{purple(level)}",
      "{msg}",
      '{gray(crayon::bold(sprintf("(%.1f MB)", pryr::mem_used() / (1024^2))))}',
      '[{gray("finished at:")} {crayon::italic(format(time, "%Y-%m-%d %H:%M:%S"))}]'
    )
  )

  logger::log_layout(local_layout, namespace = "chunk")



  logger::log_appender(
    (function(...) {
      appender <- logger::appender_file(...)
      function(lines) {
        appender(cli::ansi_strip(lines))
        logger::appender_console(lines)
      }
    })(file = "chunk_times_log.txt"),
    namespace = "chunk"
  )
}, env = environment(log_chunk_end_time))

logger::log_appender(logger::appender_tee("logger.txt"))
