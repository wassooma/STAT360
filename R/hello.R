#'stat360
#'@export
stat_start=function(NAME, alpha=0, B10=0, ..., m=0, g=0){

  data_<- read.table(NAME,header=TRUE)

  alpha<<-alpha
  #variables
  Y<<- matrix(data_$Yvar,ncol=1)

  p<<-ncol(data_)
  n<<-length(Y)

  data_<-subset(data_, select=-c(Y))
  data_<-data_[,order(names(data_))]

  X<<-matrix(c(rep(1,n)))
  for(i in 1:p-1){
    X=cbind(X,data_[,i])
  }
  # X<-X[-c(1)]
  colnames(X)=paste("X",0:(p-1),sep="")

  rm(data_)

  f1=function(){
    v_=matrix()
    for(i in 2:ncol(X)){
      v_=cbind(v_,paste0("X$",colnames(X)[i]))
    }
    return(v_[-c(1)])
  }

  fit=lm(reformulate(c(f1()),"Y"))
  print(summary(fit))

  X<<-data.matrix(X)

  J=matrix(c(rep(1,n)))
  JJ=matrix(rep(cbind(c(),J),n), ncol=n)
  I_=diag(n)

  if(nrow(matrix(c(...),ncol=1))!=(ncol(X)-1)) Xh<<-matrix(c(1,rep(0,ncol(X)-1)),ncol=1)
  else Xh<<-matrix(c(1,...),ncol=1)

  b<<-solve(t(X)%*%as.matrix(X))%*%t(X)%*%Y
  Fitted<<-as.matrix(X)%*%b

  H<<-as.matrix(X)%*%solve(t(X)%*%as.matrix(X))%*%t(X)
  # Fitted<<-H%*%Y

  #################################################

  #ANOVA
  #residuals
  residual<<-Y-Fitted

  #prediction deviations (error)
  prediction_deviation<<-Fitted-mean(Y)

  #total deviations
  total_deviation<<-Y-mean(Y)

  #sum of squares
  SSTO<<-sum(total_deviation^2)
  SSR<<-sum(prediction_deviation^2)
  SSE<<-sum(residual^2)

  #mean of squares
  MSTO<<-SSTO/(n-1)
  MSR<<-SSR/(p-1)
  MSE<<-SSE/(n-p)

  #semistudentized residuals
  e_star1<<-(residual-mean(residual))/sqrt(MSE)

  #estimated variance of residuals
  s2_residual<<-MSE*(I_-H)


  # plot(X,Y)
  # lines(X,Fitted)
  #
  # plot(sort(X),residual)
  # lines(sort(X),rep(0,n))


  #ANOVA in R

  print(anova(fit))

  #coefficient of determination (R^2) and correlation coefficient (r)
  R2<<-SSR/SSTO
  R2a<<-1-(n-1)/(n-p)*SSE/SSTO
  r<<-sqrt(R2)
  # if (b1>0) print(+r)
  # else print(-r)

  ################################################

  #estimated variance
  s2_b<<-MSE*solve(t(X)%*%as.matrix(X))
  s_b<<-sqrt(s2_b)
  # s2_b1_b2<<-(s2_b[2,3])^2

  ################################################

  statistical_test=function(letter){

    test_statistic<-0
    table_value<-0
    p_value<-0
    upper<-0
    lower<-0

    estInterval=function(){
      print(paste0("The estimated interval of ", letter, " is: (", lower, ", ", upper,")"))
    }

    conclusion=function(){
      print(paste0(letter," test results:"))

      if(test_statistic!=0 & table_value!=0){
        if(abs(test_statistic)<=table_value) print(paste0("conclude H0 because test_statistic < table_value: ", test_statistic, " < ", table_value))
        else print(paste0("conclude Ha because test_statistic > table_value: ", test_statistic, " > ", table_value))
      }

      if(test_statistic!=0 & table_value!=0){
        if(B10<=upper & B10>=lower) print(paste0("conclude H0 because B10 ",B10, " lies inside interval: (",lower , ", ", upper, ")"))
        else print(paste0("conclude Ha because B10 ",B10, " does not lie inside interval: (",lower , ", ", upper, ")"))
      }

      if(p_value!=0 & alpha!=0){
        if(p_value>alpha) print(paste0("conclude H0 because p_value is larger than significance level: ", p_value, " > ", alpha))
        else print(paste0("conclude Ha because p_value is smaller than significance level: ", p_value, " < ", alpha))
      }
    }

    #F-test for linear regression
    if(letter=="f"){
      test_statistic<-MSR/MSE
      table_value<-qf(1-alpha,p-1,n-p)
      p_value<-pf(test_statistic,p-1,n-p, lower.tail = FALSE)
      conclusion()
    }

    #t-test for Bk
    if(letter=="bk"){
      for(i in 1:ncol(s2_b)){
        bk=b[i,1]
        s2_bk=s2_b[i,i]
        s_bk=sqrt(s2_bk)
        test_statistic=(bk-B10)/s_bk
        if(B10==0){
          table_value<-qt(1-alpha/2,n-p)
          p_value<-2*pt(-abs(test_statistic),n-p)
        }
        else{
          table_value<-qt(1-alpha,n-p)
          p_value<-pt(-abs(test_statistic),n-p)
        }
        upper<-bk+table_value*s_bk
        lower<-bk-table_value*s_bk
        conclusion()
        estInterval()
      }
    }

    #t-test for Yh_hat
    if(letter=="Yh_hat"){
      test_statistic<-(Yh_hat-B10)/s_Yh_hat
      if(B10==0){
        table_value<-qt(1-alpha/2,n-p)
        p_value<-2*pt(-abs(test_statistic),n-p)
      }
      else{
        table_value<-qt(1-alpha,n-p)
        p_value<-pt(-abs(test_statistic),n-p)
      }
      upper<-Yh_hat+table_value*s_Yh_hat
      lower<-Yh_hat-table_value*s_Yh_hat
      conclusion()
      estInterval()
    }

    #F-test for lack of fit
    if(letter=="f_LF"){
      test_statistic<-MSLF/MSPE
      table_value<-qf(1-alpha, c_-p, n-c_)
      conclusion()
    }

    #bonferroni_inequality_to_estimate_Bk
    if(letter=="B_CI"){
      cat(paste0("\nBonferroni_inequality_to_estimate_Bk:\n"))
      for(i in 1:nrow(b)){
        bk=b[i,1]
        s2_bk=s2_b[i,i]
        s_bk=sqrt(s2_bk)
        table_value=qt(1-alpha/4,n-p)

        upper<-bk+table_value*s_bk;
        lower<-bk-table_value*s_bk;

        estInterval()
      }
    }
    cat("\n")
  }
  ################################################

  #Inference on Yh_hat (interval estimation)
  Yh_hat<<-(t(b)%*%Xh)[1,1]

  #estimated variance
  s2_Yh_hat<<-(MSE*t(Xh)%*%solve(t(X)%*%as.matrix(X))%*%Xh)[1,1]
  s_Yh_hat<<-sqrt(s2_Yh_hat)

  ################################################

  #F-test for linear regression
  statistical_test("f")

  #t-tests on Bk and Yh_hat
  statistical_test("bk")
  statistical_test("Yh_hat")

  ################################################

  t_table<<-qt(1-alpha/2,n-p)

  ################################################

  #Inference on pred (prediction of new observation estimation for given Xh)
  #estimated variance
  s2_pred<<-MSE+s2_Yh_hat
  s_pred<<-sqrt(s2_pred)


  #prediction interval
  upperPI_pred<<-Yh_hat+t_table*s_pred
  lowerPI_pred<<-Yh_hat-t_table*s_pred

  ################################################

  #Inference on predmean (prediction mean of m new observation for given Xh)
  #estimated variance
  s2_predmean<<-MSE/m+s2_Yh_hat
  s_predmean<<-sqrt(s2_predmean)

  #prediction interval
  upperPI_predmean<<-Yh_hat+t_table*s_predmean
  lowerPI_predmean<<-Yh_hat-t_table*s_predmean

  ################################################

  #confidence Region for the Regression Surface (PB?)
  #working-hotelling approach
  W2<<-2*qf(1-alpha,p,n-p)
  W<<-sqrt(W2)

  upperPB<<-Yh_hat+W*s_Yh_hat
  lowerPB<<-Yh_hat-W*s_Yh_hat

  ################################################

  # # general test approach
  # SSE_F_general_test<<-SSE
  # SSE_R_general_test<<-SSTO
  #
  # #F test (test statistic of {SSE(R) - SSE(F)} )
  # df_F_general_test<<-n-2
  # df_R_general_test<<-n-1
  #
  # F_star_general_test<<-((SSE_R-SSE_F)/(df_R-df_F))/((SSE_F)/df_F)
  #
  # F_table_general_test<<-qf(1-alpha, df_R-df_F, df_F)
  #
  # if(F_star_general_test<=F_table_general_test) print("conclude H0")
  # else print("conclude Ha")

  ################################################

  #exepected_value_of_e_under_normality
  k=1:n
  Ee<<-sqrt(MSE)*qnorm((k-0.375)/(n+0.25))

  plot(sort(Ee),sort(residual))

  fit2<- lm(sort(residual) ~ sort(Ee))
  anova(fit2)
  summary(fit2)

  ################################################

  #F_test_for_lack_of_fit
  data_table<<-data.frame(subset(X,select=-c(X0)),y=Y,nj=c(rep(1,n)),
                            mean=c(rep(0,n)),stDev=c(rep(NA,n)))

  for(i in 1:nrow(data_table)){
    if(data_table[i,1]==0)next
    for(j in i:nrow(data_table)){
      if(i==j) next

      Q1=subset(data_table, select=c(1:p-1))[i,]
      Q2=subset(data_table, select=c(1:p-1))[j,]
      rownames(Q1)=NULL
      rownames(Q2)=NULL

      if(identical(Q1,Q2)){
        data_table[i,"y"]<<-paste0(data_table[i,"y"],", ",data_table[j,"y"])
        data_table[i,"nj"]<<-data_table[i,"nj"]+1

        data_table[j,]<<-0
      }
    }
  }

  data_table<<-subset(data_table,X1!=0)
  rownames(data_table) <<- NULL

  c_<<-nrow(data_table)

  for(i in 1:nrow(data_table)){
    v_=c()
    for(j in 1:data_table[i,"nj"]){
      v_=cbind(v_,as.double(tokenize_string(data_table[i,"y"])
                            [1+(j-1)*3,1]))
      data_table[i,"mean"]<<-mean(v_)
      data_table[i,"stDev"]<<-sd(v_)
    }
  }

  data_table[is.na(data_table)]<<-0
  print(data_table)

  SSPE<<-0
  for(i in 1:nrow(data_table)){
    SSPE<<-SSPE+data_table[i,"stDev"]^2*(data_table[i,"nj"]-1)
  }

  MSPE<<-SSPE/(n-c_)
  SSLF<<-(SSE-SSPE)
  MSLF<<-SSLF/(c_-2)



  #lack of fit test
  f2=function(){
    z_="Y~0 + as.factor(X$X1)"
    if(ncol(X)==2) return(z_)
    for(i in 3:ncol(X)){
      z_=paste0(z_," + ","as.factor(X$",colnames(X)[i],")")
    }
    return(z_)
  }

  as.data.frame(X)#??
  Full<-lm(as.formula(f2()))
  print(anova(fit,Full))

  statistical_test("f_LF")

  ################################################

  #bonferroni_inequality_to_estimate_Bk
  B_CI<<-qt(1-alpha/4,n-p)

  statistical_test("B_CI")

  ################################################

  #working_hotteling_procedure_for_simultaneous_estimation_of_mean_response
  W2<<-p*qf(1-alpha,p,n-p)
  W<<-sqrt(W2)

  upperCI_WH_Yh_hat<<-Yh_hat+W*s_Yh_hat
  lowerCI_WH_Yh_hat<<-Yh_hat-W*s_Yh_hat

  ################################################

  #bonferroni_procedure_for_simultaneous_estimation_of_mean_response
  B_CP<<-qt(1-(alpha/(2*g)),n-p)#################it was 2*g before changing it to p*g

  upperCP_bonferroni_Yh_hat<<-Yh_hat+B_CP*s_Yh_hat
  lowerCP_bonferroni_Yh_hat<<-Yh_hat-B_CP*s_Yh_hat

  ################################################

  #scheffe_procedure_for_simultaneous_PI_for_g_new_observations
  S2<<-p*qf(1-alpha,g,n-p)#################it was 2*qf(...) before changing it to p*qf(...)
  S<<-sqrt(S2)

  upperPI_scheffe_pred<<-Yh_hat+S*s_pred
  lowerPI_scheffe_pred<<-Yh_hat-S*s_pred

  ################################################

  #bonferroni_procedure_for_simultaneous_PI_for_g_new_observations
  upperPI_bonferroni_pred<<-Yh_hat+B_CP*s_pred
  lowerPI_bonferroni_pred<<-Yh_hat-B_CP*s_pred

  ################################################

  #effect of list1 given list2
  #notation of fct is similar to inside brackets of SSR(list1|list2)
  test_multi<<-function(..., nb){

    listOfNum=matrix(c(...),ncol=1)
    list1=matrix(listOfNum[1:nb,],ncol=1)
    list2=matrix(listOfNum[(nb+1):nrow(listOfNum),],ncol=1)

    find_fct_1=function(someList){
      v_="Y ~ "
      for(i in 1:nrow(someList)){
        v_=paste0(v_, "X[,",(someList[i]+1), "]")
        if (i!=nrow(someList)) v_=paste0(v_, "+")
      }
      return(v_)
    }

    #uses the result of v_ from find_fct_2
    find_fct_2=function(someList){
      v_=paste0(find_fct_1(list2), "+")
      for(i in 1:nrow(someList)){
        v_=paste0(v_, "X[,",(someList[i]+1), "]")
        if (i!=nrow(someList)) v_=paste0(v_, "+")
      }
      return(v_)
    }

    q_<<-nrow(list2)+1

    dfF=n-p
    dfR=n-q_

    aov=anova(lm(as.formula(find_fct_1(list2))),lm(as.formula(find_fct_2(list1))))
    print(aov)

    numerator=as.matrix(aov["Sum of Sq"])[2,1]/(dfR-dfF) #effect of list1 given list2
    denominator=MSE

    F_star=numerator/denominator
    f_table=qf(1-alpha, dfR-dfF, dfF)

    p_value=pf(F_star, dfR-dfF, dfF, lower.tail = FALSE)


    cat(paste0("\ndfF and dfR are: ", dfF, ", ", dfR, " respectively\nnumerator = ",
               numerator, "\ndenominator = ", denominator, "\nF_star = ", F_star,
               "\nf_table = ", f_table, "\np-value = ", p_value))

  }

  ################################################

  r2Y<<-function(..., nb){

    listOfNum=matrix(c(...),ncol=1)
    list1=matrix(listOfNum[1:nb,],ncol=1)
    list2=matrix(listOfNum[(nb+1):nrow(listOfNum),],ncol=1)

    find_fct_1=function(someList){
      v_="Y ~ "
      for(i in 1:nrow(someList)){
        v_=paste0(v_, "X[,",(someList[i]+1), "]")
        if (i!=nrow(someList)) v_=paste0(v_, "+")
      }
      return(v_)
    }

    #uses the result of v_ from find_fct_2
    find_fct_2=function(someList){
      v_=paste0(find_fct_1(list2), "+")
      for(i in 1:nrow(someList)){
        v_=paste0(v_, "X[,",(someList[i]+1), "]")
        if (i!=nrow(someList)) v_=paste0(v_, "+")
      }
      return(v_)
    }

    aov1=anova(lm(as.formula(find_fct_1(list2))),lm(as.formula(find_fct_2(list1))))
    print(aov1)

    aov2=anova(lm(as.formula(find_fct_1(list2))))
    print(aov2)

    numerator=as.matrix(aov1["Sum of Sq"])[2,1]
    denominator=aov2["Residuals", "Sum Sq"]

    answer=numerator/denominator

    cat(paste0("\nnumerator = ", numerator, "\ndenominator = ", denominator, "\nanswer = ", answer))

  }
}


