#' mlmm function for missing response in multilevel model.
#'
#' @description mlmm function handles Bayesian multilevel model with responses that is not-missing-at-random. It is motivated from analysing mass spectrometry data with respondent dependant non-missing-at-random missingness, and when there are known covariates associating with the missingness.
#' @import MASS Matrix Rcpp stats4 ggplot2
#'
#' @param formula_completed The main regression model formula. It has the same formula format as lmr() and it is used to define the first level response and its explanatory variables.
#' @param formula_missing The logistic regression model formula. It has the same formula as formula_completed.
#' @param formula_subject The second level formula in the multilevel model which is used to define responses such as subject and its explanatory variables.
#' @param pdata The dataset contains response and predictors in a long format. Response is a vector with an indictor variable to define the corresponding unit. The data needs to have the following rudimental variables: the indicator variable for first level response, second level indicator variable for subject, or a sampling unit, an indicator for missingness and indictor of censoring. Missingness and censor are two different classification, there should not have any overlap between missingness and censored. Data structure can be referenced from the example and reference papers.
#' @param respond_dep_missing An indicator of whether response value is missing-dependant.
#' @param pidname Variable name to define the multilevel response unit , i.e. protein name or gene name
#' @param sidname Vriable name to define the subject unit, i.e. patient id or sampling id
#' @param iterno Number of iterations for the posterior samplings
#' @param chains rstan() parameter to define number of chains of posterior samplings
#' @param pathname Path to save output summary results
#' @param thin rstan() parameter to define the frequency of iterations saved
#' @param seed random seed for rstan() function
#' @param algorithm rstan() parameter which has three options c(NUTS,HMC,Fixed_param).
#' @param warmup Number of iterations for burn-out in stan.
#' @param adapt_delta_value  Adaptive delta value is an adaptation parameters for sampling algorithms,default is 0.85, value between 0-1.
#' @param savefile A logical variable to indicate if the sampling files are to be saved.
#' @param usefit A logical variable to indicate if the model use the existing fit.
#' @param stanfit The name of the fitted stan model read from .rds fle.  
#'
#' @return Return of the function is the result fitted by stan(). It will have the summarized parameters from all chains and summary results for each chain.Plot() function will return the visualization of the mean and parameters.
#' @examples
#' \dontrun{
#' library(MASS)
#' set.seed(150)
#' var2=abs(rnorm(1000,0,1));treatment=c(rep(0,500),rep(1,500))
#' geneid=rep(c(1:100),50);sid=c(rep(c(1:25),20),rep(c(26:50),20))
#' cov1=rWishart(1,df=100,Sigma=diag(rep(1,100)))
#' u=rnorm(100,0,1)
#' mu=mvrnorm(n=1,mu=u,cov1[,,1])
#' sdd=rgamma(1,shape=1,scale=1/10)
#' var1=(1/0.85)*var2+2*treatment
#' for (i in 1:1000) {var1[i]=var1[i]+rnorm(1,mu[geneid[i]],sdd)}
#' miss_logit=var2*(-0.9)+var1*(0.01)
#' probmiss=exp(miss_logit)/(exp(miss_logit)+1)
#' miss=rbinom(1000,1,probmiss);table(miss)
#' pdata=data.frame(var1,var2,treatment,miss,geneid,sid)
#' for ( i in 1:1000) if (pdata$miss[i]==1) pdata$var1[i]=NA;
#' pidname="geneid";sidname="sid";
#' #copy and paste the following formulas to the mmlm() function respectively
#' formula_completed=var1~var2+treatment
#' formula_missing=miss~var2
#' formula_censor=censor~1
#' formula_subject=~treatment
#' pathdir=getwd()
#'
#' model3=mlmm(formula_completed=var1~var2+treatment,formula_missing=miss~var2,
#' formula_subject=~treatment,pdata=pdata,respond_dep_missing=TRUE,
#' pidname="geneid",sidname="sid",pathname=pathdir,iterno=5,chains=1,savefile=FALSE)
#' }
#' @export

mlmm=function(formula_completed,formula_missing,formula_subject,pdata,respond_dep_missing=TRUE,pidname,sidname,iterno=100,chains=3,pathname,thin=1,seed=125,algorithm="NUTS",warmup=floor(iterno/2),adapt_delta_value=0.85,savefile=TRUE,usefit=T,stanfit)
{ 
  current.na.action = options('na.action')
  options(na.action='na.pass')
  t=stats::terms(formula_completed)
  mf=stats::model.frame(t,pdata,na.action='na.pass')
  mm=stats::model.matrix(mf,pdata)
  t2=stats::terms(formula_missing)
  mf2=stats::model.frame(t2,pdata,na.action='na.pass')
  mm2=stats::model.matrix(mf2,pdata)
  missing=stats::model.response(mf2)

  if (!is.null(formula_subject))
  {tt3=stats::terms(formula_subject);mf3=stats::model.frame(tt3,pdata,na.action='na.pass');mm3=stats::model.matrix(mf3,pdata)}
  y_all=stats::model.response(mf)
  options('na.action' = current.na.action$na.action)

  #######################Prepare data
  ns=length(y_all)
  nmiss=table(missing)[2]
  nobs=ns-nmiss

  datamiss=subset(pdata,(missing==1))
  dataobs=subset(pdata,(missing!=1))
  npred=dim(mm)[2]
  npred_miss=dim(mm2)[2]
  npred_sub=dim(mm3)[2]

  sid=dataobs[,sidname]
  sid_m=datamiss[,sidname]
  nsid=length(table(pdata[,sidname]))

  pid=dataobs[,pidname]
  pid_m=datamiss[,pidname]

  np=length(table(pdata[,pidname]))

  pred=mm[(missing!=1),]
  pred_miss=mm2[(missing!=1),]
  pred_sub=mm3[(missing!=1),]
  pred_m=mm[(missing==1),]
  pred_miss_m=mm2[(missing==1),]
  pred_sub_m=mm3[(missing==1),]

  miss_m=datamiss[,colnames(mf2)[1]]
  miss_obs=dataobs[,colnames(mf2)[1]]
  if (respond_dep_missing) respond_dep=1 else respond_dep=0

  #data for prior
  R=as.matrix(Matrix::Diagonal(npred))
  prec=matrix(nrow=npred,ncol=npred)
  for ( i in 1:npred)
    for (j in 1:npred)
    {if (i==j) prec[i,j]=0.1 else prec[i,j]=0.005}
  mn=rep(0,npred)
  Sigma_sd=rep(10,npred)

  prstan_data=list(nobs=nobs,nmiss=nmiss,nsid=nsid,np=np,npred=npred,npred_miss=npred_miss,npred_sub=npred_sub,
                     respond_dep=respond_dep,
                     y=y_all[(missing!=1)],sid=sid,pid=pid,pred=pred,pred_miss=pred_miss,pred_sub=pred_sub,miss_obs=miss_obs,
                     sid_m=sid_m,pid_m=pid_m,pred_m=pred_m,pred_miss_m=pred_miss_m,pred_sub_m=pred_sub_m,miss_m=miss_m,
                     R=R,Sigma_sd=Sigma_sd,prec=prec,mn=mn)
  if (respond_dep==1) parsstr=c("U","beta2","alpha","alpha_response") else parsstr=c("U","beta2","alpha")
  
  initvalue1=function () {setinitvalues(npred=npred,np=np,npred_miss=npred_miss,npred_sub=npred_sub,nmiss=nmiss,nsid=nsid)}

  if (usefit==TRUE) fitmlmm=rstan::stan(fit=stanfit,data=prstan_data,iter=iterno,init=initvalue1,pars=parsstr,seed=seed,thin=thin,algorithm=algorithm,warmup=warmup,chains=chains,control=list(adapt_delta=adapt_delta_value),sample_file=paste(pathname,"samples",sep=""))
    else            fitmlmm=rstan::stan(model_code=readstancode(modelname="mlmm"),data=prstan_data,iter=iterno,init=initvalue1,pars=parsstr,seed=seed,thin=thin,algorithm=algorithm,warmup=warmup,chains=chains,control=list(adapt_delta=adapt_delta_value),sample_file=paste(pathname,"samples",sep=""))
  print(fitmlmm)
  if (savefile==TRUE) utils::write.csv(as.array(fitmlmm),file=paste(pathname,"outsummary.csv",sep=""),row.names=TRUE)
  return(fitmlmm)
}

