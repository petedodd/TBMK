/* estimate incidence*/
functions{
  //apply log OR to probability
  real alo(real logOR, real p){
    return inv_logit(logit(p)+logOR);
  }
  matrix alom(real logOR, matrix p){//overloading only since v2.29
    return inv_logit(logit(p)+logOR);
  }
}
data{
  int N; //no countries
  int couple;//couple notifs and progression?
  real mu_styblo; real<lower=0> sig_styblo; //beta
  real mu_notes[4]; real<lower=0> sig_notes[4]; //proportion of note TBM
  real mu_prog[4]; real<lower=0> sig_prog[4]; //progrsesion risk to TBM
  real mu_age[2]; real<lower=0> sig_age[2]; //proportion of note in each age
  real<lower=0> bcgProtA; real<lower=0> bcgProtB; //BCG protection as HR
  matrix[N,2] notif_data;
  matrix[N,4] notif_dataW;
  matrix[N,4] CDR_mu;   matrix[N,4] CDR_sg; //prior for detection
  vector[N] BCG_coverage;
  matrix[N,4] POPS;
  vector[N] PREVmu; vector[N] PREVsg;
  int RN[6];//
  int whoregionAFR[RN[1]];
  int whoregionAMR[RN[2]];
  int whoregionEMR[RN[3]];
  int whoregionEUR[RN[4]];
  int whoregionSEA[RN[5]];
  int whoregionWPR[RN[6]];
  real seqA[4]; real seqB[4]; //meta-analysis sequelae
  real CFRtxA[4]; real CFRtxB[4]; //meta-analysis CFR
  real lgOR[4]; real lgORsd[4];  //CFR age ORs from ZAF data
  real hiv_lmn[N]; real hiv_lsg[N]; //HIV in TB prevalence
  real mu_hivd; real sig_hivd; //OR for mortality given HIV
  real mu_hivi; real sig_hivi; //OR for TBM given HIV
}
transformed data{

}
parameters{
  vector<lower=0>[N] PREV;
  real<lower=0> styblo;
  real<lower=0,upper=1> BCGprot;
  //matrix<lower=0.0001,upper=0.999>[N,4] propTBM;
  vector<lower=0.0001,upper=0.999>[4] propTBM;
  //vector<lower=0,upper=1>[4] propTBM;
  vector<lower=0.01,upper=0.999>[2] propAge;
  vector<lower=0,upper=1>[4] prog;
  real lgcfrOR[4];
  real hivdOR;//log OR
  real hiviOR;//log OR
  vector<lower=0,upper=1>[N] hivintb;
}
transformed parameters{
  vector<lower=0>[N] foi;
  vector<lower=0>[N] tbivec;
  matrix<lower=0>[N,4] tbimat;
  vector<lower=0,upper=1>[4] propAgeLong;
  matrix<lower=0,upper=1>[N,4] propAmat;
  matrix[N,4] propTBMmat;
  matrix[N,4] propTBMmat0;
  matrix[N,4] propTBMmat1;
  matrix<lower=0>[N,4] TBMnotes;
  matrix<lower=0>[N,4] TBMI;
  matrix<lower=0,upper=1>[N,4] progmat;
  matrix<lower=0>[N,4] NoI;
  matrix[N,4] hivmat;//HIV matrix
  hivmat = rep_matrix(hivintb,4);
  // //force of infection
  foi = PREV * styblo * 0.5/1e5;
  propTBMmat0 = rep_matrix(propTBM',N);// .* hivnmat;// HIV -ve
  propTBMmat1 = alom(hiviOR,propTBMmat0);// .* hivmat;
  propTBMmat = propTBMmat0 .* (1.0 - hivmat) + propTBMmat1 .* hivmat;
  propAgeLong[1] = propAge[1];
  propAgeLong[3] = propAge[2];
  propAgeLong[2] = 1-propAge[1];
  propAgeLong[4] = 1-propAge[2];
  propAmat =  rep_matrix(propAgeLong',N);
  TBMnotes = propAmat .* propTBMmat .* notif_dataW;//
  tbivec = foi .* (1-(1-BCGprot) * 1e-2 * BCG_coverage);//NOTE
  tbimat = rep_matrix(tbivec,4);
  progmat = rep_matrix(prog',N);
  TBMI = tbimat .* progmat .* POPS;
  // CDR estimate
  NoI = TBMnotes ./ TBMI;
}
model{
  PREV ~ lognormal(PREVmu,PREVsg);
  hivintb ~ lognormal(hiv_lmn,hiv_lsg);
  styblo ~ lognormal(mu_styblo,sig_styblo);
  BCGprot ~ beta(bcgProtA,bcgProtB);
  propTBM ~ lognormal(mu_notes,sig_notes);//
  propAge ~ lognormal(mu_age,sig_age);//also bad but order less //(-1,0.1)
  //propAge ~ beta(mu_age,sig_age);//NOTE debug
  prog ~ lognormal(mu_prog,sig_prog);
  //ORs
  lgcfrOR ~ normal(lgOR,lgORsd); //CFRs age
  hivdOR ~  normal(mu_hivd, sig_hivd); //log OR for mortality given HIV
  hiviOR ~ normal(mu_hivi, sig_hivi); //log OR for TBM given HIV
  if(couple>0){
      target += lognormal_lpdf( to_vector(NoI) | to_vector(CDR_mu), to_vector(CDR_sg) );
  }
}
generated quantities{
  /* age & region specific HIV proportions */
  matrix<lower=0,upper=1>[6,4] regional_fracHIVinc; matrix<lower=0,upper=1>[6,4] regional_fracHIVmort;
  vector<lower=0,upper=1>[4] global_fracHIVincAge; vector<lower=0,upper=1>[4] global_fracHIVmortAge;
  matrix<lower=0,upper=1>[N,4] hivinTBM =  propTBMmat1 .* hivmat ./ propTBMmat;//HIV prev in TBM
 //reformed CFRs
  real CFRtxR[4] = beta_rng(CFRtxA,CFRtxB); //CFR treated
  real CFRtxRadj[4] = {alo(lgcfrOR[1],CFRtxR[1]),
    alo(lgcfrOR[2],CFRtxR[2]),
    alo(lgcfrOR[3],CFRtxR[3]),
    alo(lgcfrOR[4],CFRtxR[4])}; //adjusted by ORs
  row_vector[4] CFRtxV = to_row_vector(CFRtxRadj);
  matrix[N,4] CFRtx = rep_matrix(CFRtxV,N); //matrix tx CFR
  matrix[N,4] CFRtxH = alom(hivdOR,CFRtx); //matrix tx CFR for HIV
  //reformed CMRs - morbidity
  real CMRtxR[4] = beta_rng(seqA,seqB);
  real CMRtxRadj[4] = {alo(lgcfrOR[1],CMRtxR[1]),
    alo(lgcfrOR[2],CMRtxR[2]),
    alo(lgcfrOR[3],CMRtxR[3]),
    alo(lgcfrOR[4],CMRtxR[4])}; //adjusted by ORs
  row_vector[4] CMRtxV = to_row_vector(CMRtxRadj);
  matrix[N,4] CMRtx = rep_matrix(CMRtxV,N);
  //calculations
  matrix<lower=0>[N,4] untreated = fdim(TBMI,TBMnotes); //caps over dx
  matrix<lower=0>[N,4] deaths = untreated + TBMnotes .* (CFRtx .* (1.0-hivinTBM) + CFRtxH .* hivinTBM);
  matrix<lower=0>[N,4] hivinc = TBMI .* hivinTBM;
  matrix<lower=0>[N,4] hivdeaths = untreated .* hivinTBM + TBMnotes .* CFRtxH .* hivinTBM;
  matrix<lower=0>[N,4] morbs = CMRtx .* (1-CFRtx .* (1.0-hivinTBM) + CFRtxH .* hivinTBM) .* TBMnotes;
  real global_mnotes = sum(TBMnotes);
  real global_minc = sum(TBMI);
  real global_deaths = sum(deaths);
  real global_morbs = sum(morbs);
  real global_fracDeathsUntreated = sum(untreated)/global_deaths ;
  real global_fracHIVInc = sum(hivinc)/global_minc;
  real global_fracHIVDeaths = sum(hivdeaths)/global_deaths ;
  
  //notes
  real regional_AFR_mnotes = sum(TBMnotes[whoregionAFR,1:4]);
  real regional_AMR_mnotes = sum(TBMnotes[whoregionAMR,1:4]);
  real regional_EMR_mnotes = sum(TBMnotes[whoregionEMR,1:4]);
  real regional_EUR_mnotes = sum(TBMnotes[whoregionEUR,1:4]);
  real regional_SEA_mnotes = sum(TBMnotes[whoregionSEA,1:4]);
  real regional_WPR_mnotes = sum(TBMnotes[whoregionWPR,1:4]);
  //incidence
  real regional_AFR_minc = sum(TBMI[whoregionAFR,1:4]);
  real regional_AMR_minc = sum(TBMI[whoregionAMR,1:4]);
  real regional_EMR_minc = sum(TBMI[whoregionEMR,1:4]);
  real regional_EUR_minc = sum(TBMI[whoregionEUR,1:4]);
  real regional_SEA_minc = sum(TBMI[whoregionSEA,1:4]);
  real regional_WPR_minc = sum(TBMI[whoregionWPR,1:4]);
  //deaths
  real regional_AFR_deaths = sum(deaths[whoregionAFR,1:4]);
  real regional_AMR_deaths = sum(deaths[whoregionAMR,1:4]);
  real regional_EMR_deaths = sum(deaths[whoregionEMR,1:4]);
  real regional_EUR_deaths = sum(deaths[whoregionEUR,1:4]);
  real regional_SEA_deaths = sum(deaths[whoregionSEA,1:4]);
  real regional_WPR_deaths = sum(deaths[whoregionWPR,1:4]);
  //morbs
  real regional_AFR_morbs = sum(morbs[whoregionAFR,1:4]);
  real regional_AMR_morbs = sum(morbs[whoregionAMR,1:4]);
  real regional_EMR_morbs = sum(morbs[whoregionEMR,1:4]);
  real regional_EUR_morbs = sum(morbs[whoregionEUR,1:4]);
  real regional_SEA_morbs = sum(morbs[whoregionSEA,1:4]);
  real regional_WPR_morbs = sum(morbs[whoregionWPR,1:4]);
  //incidence - by age
  real age_u1_minc = sum(TBMI[,1]);
  real age_1to4_minc = sum(TBMI[,2]);
  real age_5to9_minc = sum(TBMI[,3]);
  real age_10to14_minc = sum(TBMI[,4]);
  //deaths - by age
  real age_u1_deaths = sum(deaths[,1]);
  real age_1to4_deaths = sum(deaths[,2]);
  real age_5to9_deaths = sum(deaths[,3]);
  real age_10to14_deaths = sum(deaths[,4]);
  //morbs - by age
  real age_u1_morbs = sum(morbs[,1]);
  real age_1to4_morbs = sum(morbs[,2]);
  real age_5to9_morbs = sum(morbs[,3]);
  real age_10to14_morbs = sum(morbs[,4]);
  /* region x age fractions HIV */
  for(i in 1:4){
    /* incidence/morbidity */
    regional_fracHIVinc[1,i] = sum(hivinc[whoregionAFR,i]) / sum(TBMI[whoregionAFR,i]);
    regional_fracHIVinc[2,i] = sum(hivinc[whoregionAMR,i]) / sum(TBMI[whoregionAMR,i]);
    regional_fracHIVinc[3,i] = sum(hivinc[whoregionEMR,i]) / sum(TBMI[whoregionEMR,i]);
    regional_fracHIVinc[4,i] = sum(hivinc[whoregionEUR,i]) / sum(TBMI[whoregionEUR,i]);
    regional_fracHIVinc[5,i] = sum(hivinc[whoregionSEA,i]) / sum(TBMI[whoregionSEA,i]);
    regional_fracHIVinc[6,i] = sum(hivinc[whoregionWPR,i]) / sum(TBMI[whoregionWPR,i]);
    /* mortality */
    regional_fracHIVmort[1,i] = sum(hivdeaths[whoregionAFR,i]) / sum(deaths[whoregionAFR,i]);
    regional_fracHIVmort[2,i] = sum(hivdeaths[whoregionAMR,i]) / sum(deaths[whoregionAMR,i]);
    regional_fracHIVmort[3,i] = sum(hivdeaths[whoregionEMR,i]) / sum(deaths[whoregionEMR,i]);
    regional_fracHIVmort[4,i] = sum(hivdeaths[whoregionEUR,i]) / sum(deaths[whoregionEUR,i]);
    regional_fracHIVmort[5,i] = sum(hivdeaths[whoregionSEA,i]) / sum(deaths[whoregionSEA,i]);
    regional_fracHIVmort[6,i] = sum(hivdeaths[whoregionWPR,i]) / sum(deaths[whoregionWPR,i]);
  }
  /* global  */
  /* inc */
  global_fracHIVincAge[1] = sum(hivinc[,1]) / age_u1_minc;
  global_fracHIVincAge[2] = sum(hivinc[,2]) / age_1to4_minc;
  global_fracHIVincAge[3] = sum(hivinc[,3]) / age_5to9_minc;
  global_fracHIVincAge[4] = sum(hivinc[,4]) / age_10to14_minc;
  /* mort */
  global_fracHIVmortAge[1] = sum(hivdeaths[,1]) / age_u1_deaths;
  global_fracHIVmortAge[2] = sum(hivdeaths[,2]) / age_1to4_deaths;
  global_fracHIVmortAge[3] = sum(hivdeaths[,3]) / age_5to9_deaths;
  global_fracHIVmortAge[4] = sum(hivdeaths[,4]) / age_10to14_deaths;
}
