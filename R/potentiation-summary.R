#' Conclusion text of a (de)potentiation analysis
#'
#' Turns the comparisons of [analyse_potentiation()] into one or two sentences.
#' Series are taken as increasing stressor-2 level in file order. An onset
#' ("... from Level 2 on") is named only when the significant series are exactly
#' a trailing run of series that all shift in the same direction; otherwise the
#' significant series are listed per direction and the text says the effect does
#' not increase consistently with the level. A shift whose EC50 lies outside the
#' tested concentrations is marked as extrapolated; so is a reference EC50.
#' @param res Result of [analyse_potentiation()] (needs `reference`,
#'   `comparisons` and optionally `curves`).
#' @param stressor1 Display name of stressor 1.
#' @return A single character string.
#' @keywords internal
#' @noRd
potentiation_summary <- function(res, stressor1 = "Stressor 1") {
  cmp <- res$comparisons
  S <- stressor1
  item <- function(i)
    paste0(cmp$series[[i]], " (EC50 ", format_num(cmp$ec50_reference[[i]]),
           " → ", format_num(cmp$ec50_series[[i]]), ", ",
           format_p(cmp$p[[i]]),
           if (isTRUE(cmp$ec50_extrapolated[[i]])) ", extrapolated", ")")
  join <- function(x) if (length(x) < 2) x else
    paste(paste(x[-length(x)], collapse = ", "), "and", x[[length(x)]])
  sig <- which(cmp$significant)
  m <- nrow(cmp)

  text <- if (m == 1) {
    if (length(sig))
      paste0(S, " is significantly ", cmp$direction[[1]], " with ", item(1), ".")
    else
      paste0("No significant influence of ", cmp$series[[1]],
             " on the toxicity of ", S, " (", format_p(cmp$p[[1]]), ").")
  } else if (!length(sig)) {
    paste0("No significant influence on the toxicity of ", S, " at any level.")
  } else {
    dirs <- unique(cmp$direction[sig])
    suffix <- identical(sig, seq(sig[[1]], m))
    if (suffix && length(dirs) == 1) {
      paste0(S, " is significantly ", dirs, " from ",
             cmp$series[[sig[[1]]]], " on: ",
             paste(vapply(sig, item, character(1)), collapse = "; "), ".")
    } else {
      parts <- vapply(intersect(c("more toxic", "less toxic"), dirs),
        function(d) paste0(S, " is significantly ", d, " with ",
                           join(vapply(sig[cmp$direction[sig] == d], item,
                                       character(1)))),
        character(1))
      text <- paste(parts, collapse = "; ")
      if (length(sig) < m)
        text <- paste0(text, "; no significant influence with ",
                       paste(cmp$series[-sig], collapse = ", "))
      paste0(text, ". The effect does not increase consistently with the level.")
    }
  }

  if (isTRUE(res$curves[[res$reference]]$ec50_extrapolated))
    text <- paste0(text, " The EC50 of ", res$reference,
                   " lies outside its tested concentrations, so every shift ",
                   "is extrapolated.")
  text
}

format_p <- function(p)
  if (p < 1e-4) "p < 0.0001" else
    paste("p =", format(signif(p, 2), scientific = FALSE))

format_num <- function(x) format(signif(x, 3), scientific = FALSE)
