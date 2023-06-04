#' envMeanPara_plotter
#'
#' @param data pop_cor file
#' @param dap_x depend on your file
#' @param dap_y depend on your file
#' @param p defualt p=1
#'
#' @return combined pop_cor plot
#' @export
#'
#' @examples Exhaustive_plotter(input=pop_cor,paras="PTS")
Exhaustive_plotter<-function(data,dap_x, dap_y,p){
pop_cors<-as.data.frame(data)
pop_cor <- subset(pop_cors, pop_cors$pop_code == p);
layout(matrix(c(1:(2*nParas)), 2, nParas, byrow = T))
for (k in 1:(2*nParas)) {
  pop_cor_0 <- subset(pop_cors, pop_cors$pop_code == p);
  pop_cor <- pop_cor_0[,c(1:4, k + 4)];
  colnames(pop_cor)[5] <- 'R';
  xs <- pop_cor$Day_x;
  ys <-  pop_cor$Day_y;
  mid_R <- median(pop_cor$R);
  pop_cor <- pop_cor[order(pop_cor$R),];
  cell_col <- floor(pop_cor$R * 12) + 13;
  pop_cor$cell_col <- cell_col;
  pop_cor_6 <- subset(pop_cor, pop_cor$window > 6); max_R <- pop_cor_6[which.max(pop_cor_6$R)[1], ];
  par(mar = c(0.5, 1.0, 1, 0.5) , mgp = c(0.05, 0.1, 0), tck = -0.01, bty = "n");
  plot(-50, -50, xlim = c(0, dap_x), ylim = c(0, dap_y), col = "white",  xlab = '', xaxt = "n", yaxt = "n", ylab = 'Days after planting', bty = "n", cex.lab = .4);
  arrows(-1, 10, -1, dap_y - 10, length = 0.05, angle = 15, lwd = .5,  col = "grey59");
  mtext(c(1, 50, 100, dap_y), side = 2, at = c(1,50, 100, dap_y), line = -1, cex = .6)

  rect(xs - 0.5, ys - 0.5, xs + 0.5, ys + 0.5, col = col_palette[pop_cor$cell_col], border = "NA")
  rect(max(pop_cor$Day_x) - 0.5, max(pop_cor$Day_y) - 0.5, max(pop_cor$Day_x) + 0.5, max(pop_cor$Day_y) + 0.5, border = "NA", col = "white", lwd = 0.001)
  #   legend("bottom", Div_Fnd_lab, bty = "n", cex = .6);

  arrows(10, dap_y + 4, dap_x - 10, dap_y + 4, angle = 15, length = 0.05, lwd = .5, col = "grey59")
  mtext("Days after planting", side = 3, at = dap_x / 2, line = -0.1, cex = .4)
  mtext(c(1, 50, 100, dap_y), side = 3, at = c(1, 50, 100, dap_y), line = -1.1, cex = .6)
  arrows(max_R$Day_x + 4,  max_R$Day_y - 4,  max_R$Day_x,  max_R$Day_y, length = 0.05, angle = 15, lwd = .5, col = "grey59")

  box_ys <- seq(1, 50, by = 2); box_xs <- rep(dap_x - 15, 25);
  rect(box_xs - .5 * 2, box_ys - 0.5 * 2, box_xs + 0.5 * 2, box_ys + 0.5 * 2, border = "NA", col = col_palette)
  text(dap_x - 10 - 5, 52, 'r', cex = .5);
  r_lab_top <- 1; r_lab_mid <- 0; r_lab_bottom <- -1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", max_R$R), sep = '');
  if (k > nParas) { r_lab_top <- -1; r_lab_bottom <- 1; max_r_lab <- paste( 'r = ', sprintf( "%.3f", 0 - max_R$R), sep = ''); mtext(side = 1, Paras[k - nParas ], line= -0.5,  cex = .5, bty = "n")}
  legend(max_R$Day_x - 4 , max_R$Day_y - 4 , c(paste( max_R$Day_x, ' to ', max_R$Day_y, ' DAP', sep = ''), max_r_lab),  cex = .6, bty = "n");
  text(dap_x - 10 + 3, 50, r_lab_top, cex = .5)
  text(dap_x - 10 + 3, 27, r_lab_mid, cex = .5);
  text(dap_x - 10 + 3, 1,  r_lab_bottom, cex = .5)
  if (length(trait$FTgdd) > 1) {
    boxplot(trait$FTgdd,   at = 145,  add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
    boxplot(trait$FTgdd,   at = 1, horizontal = T, add = TRUE, xaxt = "n", yaxt = "n", xlab = '', ylab = '', width = 10, pch = 19, cex = .3, boxwex = 4, lwd = .4, col = "gold", border = "grey");
    text(mean(trait$FTgdd), 5, 'Days to anthesis', cex = .5)
    text(mean(trait$FTgdd), 10, paste('Trait: ', trait, sep = ''), cex = .6)
  }
}
}
