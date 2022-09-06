get_lum <- function(col) {
  #https://stackoverflow.com/questions/596216/formula-to-determine-perceived-brightness-of-rgb-color

  # list vs vector?!
  if (length(col) == 1 && methods::is(col, "character")) {
    col <- unname(grDevices::col2rgb(col)[,1])
  }
  if (!is.numeric(col) || length(col) != 3) {
    stop("col has to be one character (hex code or name of a color), or a numeric of length 3 (rgb).")
  }

  col = col / 255
  col = .RGBtoLin(col)
  Y = (0.2126 * col[1] + 0.7152 * col[2] + 0.0722 * col[3])
  Ls <- .YtoLstar(Y)
  return(Ls)
}

.RGBtoLin <- function(cc) {
  sapply(cc, function(x) {
    if (x <= 0.04045) {
      x / 12.92
    } else {
      ((x + 0.055)/1.055)^2.4
    }
  })
}

.YtoLstar <- function(Y) {
  sapply(Y, function(x) {
    if (x <= (216/24389)) {
      x * (24389/27)
    } else {
      x^(1/3) * 116 - 16
    }
  })
}
