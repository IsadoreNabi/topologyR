# Comprueba que cada entrada de los bloques @references de R/*.R coincide, salvo el corte de renglones, con una
# entrada de la lista canónica dev/REFERENCIAS_APA7.md, y que toda entrada canónica se usa al menos una vez.
# Uso, desde la raíz del paquete: Rscript dev/verificar_referencias.R   (código de salida 0 si no hay defectos)
squash <- function(x) gsub("\\s+", " ", trimws(paste(x, collapse = " ")))
md <- readLines("dev/REFERENCIAS_APA7.md", encoding = "UTF-8")
canon <- character(0); cur <- character(0); inblock <- FALSE
for (ln in c(md, "")) {
  if (grepl("^    \\S", ln) || (inblock && grepl("^      \\S", ln))) {
    cur <- c(cur, ln); inblock <- TRUE
  } else {
    if (length(cur)) canon <- c(canon, squash(cur))
    cur <- character(0); inblock <- FALSE
  }
}
usadas <- setNames(integer(length(canon)), canon)
defectos <- 0L; total <- 0L
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  L <- readLines(f, encoding = "UTF-8")
  starts <- grep("^#' @references\\s*$", L)
  for (s in starts) {
    j <- s + 1L; block <- character(0)
    while (j <= length(L) && grepl("^#'", L[j]) && !grepl("^#' @", L[j])) { block <- c(block, sub("^#' ?", "", L[j])); j <- j + 1L }
    entradas <- split(block, cumsum(!nzchar(trimws(block))))
    for (e in entradas) {
      e <- e[nzchar(trimws(e))]
      if (!length(e)) next
      total <- total + 1L
      k <- squash(e)
      if (k %in% canon) usadas[k] <- usadas[k] + 1L else {
        defectos <- defectos + 1L
        cat("NO CANONICA en ", f, ", renglón ", s, ":\n  ", k, "\n", sep = "")
      }
    }
  }
}
sin_uso <- names(usadas)[usadas == 0L]
for (k in sin_uso) cat("CANONICA SIN USO: ", k, "\n", sep = "")
cat(sprintf("entradas en bloques: %d | canónicas: %d | no canónicas: %d | canónicas sin uso: %d\n",
            total, length(canon), defectos, length(sin_uso)))
quit(status = if (defectos == 0L && length(sin_uso) == 0L) 0L else 1L)
