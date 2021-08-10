within_variance_comp_emmeans = function(ind,peptide,data,metacond){
  tempdata <- cbind(reponse=data[peptide,,ind],data.frame(fact=factor(metacond)))
  templm <- lm(reponse~fact,data=tempdata)
  return(vcov((emmeans::emmeans(templm,specs = "fact"))))
}