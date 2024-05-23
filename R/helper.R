dt_table <- function(data.frame, title="", filename="", rownames=F){
  
  table <- DT::datatable(data.frame, options = list(dom = 'Bfrtip', scrollX=T,
                                                    buttons =
                                                      list(list(extend='copy',title=filename),
                                                           list(extend='csv',      title=filename),
                                                           list(extend='excel', filename=filename),
                                                           list(extend='pdf',   filename=filename),
                                                           list(extend='print', filename=filename)),
                                                    autoWidth = F, pageLength = 5))
  return(table)
}