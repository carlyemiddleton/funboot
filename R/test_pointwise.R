#This function performs the pointwise test for colocalization at a specific radius

test_pointwise <- function(data=data, from.cell, to.cell, qc.cellcount.cutoff, P, perm.yn=F, R=200, inc=1,
                           image.dims, summary.function='L'){

  #data is any data.frame with the following named columns: mage_number, cell_metacluster, cell_x, cell_y
  #image_number must be unique, not patient1 image1, patient2 image1, etc.  also best to make cell_id unique

}
