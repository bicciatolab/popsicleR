.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0(cat("Hello and welcome to PoPsicle, an interactive pipeline for the preprocessing of single cell data.\nIn this workflow, messages are colour coded: \n"),
           cat(green("green messages will provide information on the ongoing step \n")),
           cat(cyan("cyan messages will require the user to provide an input for advancing the script \n")),
           cat(red("red messages will report missing information or wrong inputs \n")),
           cat(silver("grey messages report functions and software system outputs. \n")))
  )
}
