# Renders a validation_report() data frame as a standalone HTML page. The point
# of the page is that it is generated, not written: every figure in it comes
# from the run that produced it, and every workbook figure comes from the same
# oracle table the test suite asserts against. Kept free of external assets so
# the file can be opened, mailed, or served from Pages as-is.

.vh_esc <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Numbers are the content here, so they are formatted for comparison rather than
# for compactness: enough significant digits to see a disagreement, no more.
.vh_num <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  if (is.na(n)) return(.vh_esc(x))
  a <- abs(n)
  trimws(if (a >= 1e6) formatC(n, format = "d", big.mark = "\u202f")
         else if (a >= 100) formatC(n, format = "f", digits = 2, big.mark = "\u202f")
         else if (a >= 1) formatC(n, format = "f", digits = 4)
         else formatC(n, format = "g", digits = 5))
}

.vh_kind_label <- c(
  match        = "reproduces workbook",
  divergent    = "deliberate difference",
  pinned       = "fixed in both",
  unidentified = "not identified by the data",
  bounded      = "inequality only")

# The headroom bar: how much of the allowed tolerance a match actually used.
# A quantity that lands at 3% of its tolerance is a far stronger result than one
# that scrapes in at 98%, and a table of numbers alone hides that difference.
.vh_bar <- function(rel, tol) {
  if (is.na(rel) || is.na(tol) || tol <= 0)
    return('<div class="bar bar--na" aria-hidden="true"></div>')
  frac <- min(rel / tol, 1)
  pct  <- max(frac * 100, 1.2)
  over <- rel > tol
  sprintf(paste0('<div class="bar%s"><span class="bar__fill" style="width:%.1f%%">',
                 '</span></div><span class="bar__pct">%s</span>'),
          if (over) " bar--over" else "", pct,
          if (over) "over" else sprintf("%.0f%% of tol.", frac * 100))
}

#' Render a validation report as a standalone HTML page
#'
#' @param report A [validation_report()] data frame.
#' @param path File to write. Defaults to `validation-report.html` in the
#'   working directory.
#' @param generated Timestamp shown in the header.
#' @param standalone Wrap the page in a full HTML document (the default). Set
#'   `FALSE` to emit just the title, styles and body content, for embedding in a
#'   host that supplies its own document scaffold.
#' @return `path`, invisibly.
#' @export
validation_report_html <- function(report, path = "validation-report.html",
                                   generated = Sys.time(), standalone = TRUE) {
  reg <- validation_datasets()
  n_match <- sum(report$kind == "match")
  n_pass  <- sum(report$kind == "match" & report$verdict == "pass")
  n_fail  <- sum(report$verdict == "fail")
  n_div   <- sum(report$kind == "divergent")

  panels <- vapply(unique(report$dataset), function(ds) {
    d <- report[report$dataset == ds, , drop = FALSE]
    title <- reg[[ds]]$title %||% ds
    rows <- vapply(seq_len(nrow(d)), function(i) {
      wb  <- if (nzchar(as.character(d$workbook[i])) && !is.na(d$workbook[i]))
               .vh_num(d$workbook[i]) else "\u2014"
      obs <- if (!is.na(d$observed[i])) .vh_num(d$observed[i]) else "not produced"
      sprintf(paste0(
        '<div class="row row--%s">',
        '<div class="row__q">%s<span class="row__param">%s</span></div>',
        '<div class="row__v"><span class="row__cap">workbook</span>%s</div>',
        '<div class="row__v"><span class="row__cap">mixdra</span>%s</div>',
        '<div class="row__d">%s</div>',
        '<div class="row__k">%s</div>',
        '</div>%s'),
        .vh_esc(d$kind[i]), .vh_esc(d$label[i]), .vh_esc(d$parameter[i]),
        wb, obs, .vh_bar(d$rel_diff[i], suppressWarnings(as.numeric(d$tolerance[i]))),
        .vh_esc(.vh_kind_label[[d$kind[i]]] %||% d$kind[i]),
        if (nzchar(d$note[i]))
          sprintf('<p class="note"><span>%s</span>%s</p>',
                  .vh_esc(d$source[i]), .vh_esc(d$note[i]))
        else if (nzchar(d$source[i]))
          sprintf('<p class="note"><span>%s</span></p>', .vh_esc(d$source[i]))
        else "")
    }, character(1))
    sprintf('<section class="panel"><h2>%s</h2><p class="panel__id">%s</p>%s</section>',
            .vh_esc(title), .vh_esc(ds), paste(rows, collapse = "\n"))
  }, character(1))

  # Placeholder substitution rather than sprintf: the template is mostly CSS,
  # and CSS is full of literal percent signs.
  html <- .vh_template()
  subs <- list(c("{{GENERATED}}", format(generated, "%d %B %Y, %H:%M")),
               c("{{PASS}}",   as.character(n_pass)),
               c("{{MATCH}}",  as.character(n_match)),
               c("{{DIV}}",    as.character(n_div)),
               c("{{FAIL}}",   as.character(n_fail)),
               c("{{PANELS}}", paste(panels, collapse = "\n")))
  for (kv in subs) html <- sub(kv[[1]], kv[[2]], html, fixed = TRUE)

  # The template is a fragment: title, styles and content, with a marker where
  # the document scaffold would close its head. A host that supplies its own
  # scaffold (an artifact viewer, a vignette) takes it as-is; a file on disk
  # needs the scaffold wrapped around it.
  html <- sub("<!--BODY-->", if (standalone) "</head>\n<body>" else "",
              html, fixed = TRUE)
  if (standalone)
    html <- paste0('<!doctype html>\n<html lang="en">\n<head>\n',
                   '<meta charset="utf-8">\n',
                   '<meta name="viewport" content="width=device-width, ',
                   'initial-scale=1">\n', html, "\n</body>\n</html>\n")
  writeLines(html, path, useBytes = TRUE)
  invisible(path)
}
