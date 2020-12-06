
#'stat360
#'@export
stat_start=function(NAME, num1=0, num2=0, num3=0, num3b=0, num4=0, num5=0, multivariable=FALSE, matrix=FALSE){
  data1<- read.table(NAME,header=TRUE)

  if(matrix==TRUE){
    #variables

    X1<<-data1$X1var
    Yvar<<- data1$Yvar
    n<<-length(Yvar)
    p<<-2
    X2<<-matrix(c(rep(0,n)),ncol=1)

    if(multivariable==TRUE){
      X2<<- data1$X2var
    }

    Y<<-matrix(Yvar, ncol=1)
    X<<-matrix(c(rep(1,n)))
    X<<-cbind(X,X1)
    if(multivariable==TRUE)X<<-cbind(X,X2)

    fit<<-lm(Y~X1+X2)
    print(summary(fit))

    J<<-matrix(c(rep(1,n)))
    JJ<<-matrix(rep(cbind(c(),J),n), ncol=n)
    I_<<-diag(n)

    Xh<<-matrix(c(1,num3),ncol=1)
    if(multivariable==TRUE)Xh<<-matrix(c(1,num3,num3b),ncol=1)

    #variables assignment
    alpha<<-num1
    B10<<-num2
    Xh1<<-num3
    Xh2<<-num3b
    m<<-num4
    g<<-num5

    #calculations
    # Sx1y<<-sum(X1*Y)-(sum(X1)*sum(Y))/n
    # Sx2x<<-sum(X2*Y)-(sum(X2)*sum(Y))/n
    # Sx1x2<<-sum(X1*X2)-(sum(X1)*sum(X2))/n
    #

    #regression line coefficients

    b<<-solve(t(X)%*%X)%*%t(X)%*%Y
    b0<<-unname(b[1,1])
    b1<<-unname(b[2,1])
    if(multivariable==TRUE)b2<<-unname(b[3,1])
    Fitted<<-b0+b1*X1
    if(multivariable==TRUE)Fitted<<-b0+b1*X1+b2*X2

    #residuals
    H<<-X%*%solve(t(X)%*%(X))%*%t(X)
    Fitted<<-H%*%Y
    residual<<-Y-Fitted

    SSE<<-sum(residual^2)
    MSE<<-SSE/(n-2)

    #########
    s2_residual<<-MSE*(I_-H)
    ###########

    # plot(X,Y)
    # lines(X,Fitted)
    #
    # plot(sort(X),residual)
    # lines(sort(X),rep(0,n))

    #estimation of residual variance sigma_Sq
    s2<<- MSE
    s<<- sqrt(MSE)

    #################################################

    #ANOVA

    #semistudentized residuals
    e_star1<<-(residual-mean(residual))/sqrt(MSE)

    #prediction deviation (error)
    prediction_deviation<<-Fitted-mean(Y)

    #total deviation
    total_deviation<<-Y-mean(Y)

    #sum of squares
    SSTO<<-sum(total_deviation^2)
    SSR<<-sum(prediction_deviation^2)
    SSE<<-sum(residual^2)

    #mean of squares
    MSTO<<-SSTO/(n-1)
    MSR<<-SSR/1
    MSE<<-SSE/(n-2)

    #F-test for linear regression
    F_star<<-MSR/MSE
    F_table<<-qf(1-alpha,p-1,n-p)
    if(F_star<=F_table) print("conclude H0")
    else print("conclude Ha")

    #t-test (using F) for linear regression
    t_star_using_F<<-sqrt(F_star)
    t_table_using_F<<-qt(1-alpha/2,n-2)
    if(t_star_using_F<=t_table_using_F) print("conclude H0")
    else print("conclude Ha")

    #ANOVA in R

    print(anova(fit))

    #coefficient of determination (R^2) and correlation coefficient (r)
    R_2<<-SSR/SSTO
    r<<-sqrt(R_2)
    if (b1>0) print(+r)
    else print(-r)

    ################################################

    #inference on B1

    #estimated variance
    s2_b<<-MSE*solve(t(X)%*%X)

    s2_b0<<-s2_b[1,1]
    s_b0<<-sqrt(s2_b0)

    s2_b1<<-s2_b[2,2]
    s_b1<<-sqrt(s2_b1)
    if(multivariable==TRUE)s2_b2<<-s2_b[3,3]
    if(multivariable==TRUE)s2_b1_b2<<-(s2_b[2,3])^2
    if(multivariable==TRUE)s_b2<<-sqrt(s2_b2)

    #t-test concerning B1
    if(B10==0){
      t_star<<-(b1-B10)/s_b1
      t_table<<-qt(1-alpha/2,n-p)

      p_value<<-2*pt(-abs(t_star),n-p)

    }
    else{
      t_star<<-(b1-B10)/s_b1
      t_table<<-qt(1-alpha,n-p)

      p_value<<-pt(-abs(t_star),n-p)
    }

    #confidence interval of B1
    upperCI_b1<<-b1+t_table*s_b1
    lowerCI_b1<<-b1-t_table*s_b1

    if(abs(t_star)<=t_table) print("conclude H0")
    else print("conclude Ha")

    if(B10<=upperCI_b1 & B10>=lowerCI_b1) print("conclude H0 because B10 lies inside CI")
    else print("conclude Ha because B10 does not lie inside CI")

    if(p_value>alpha) print("conclude H0 because p_value is larger than significance level: alpha")
    else print ("conclude Ha because p_value is smaller than significance level: alpha")

    t_star_b1<<-t_star
    t_table_b1<<-t_table
    p_value_b1<<-p_value

    ################################################

    #inference on B0

    #confidence interval
    upperCI_b0<<-b0+t_table*s_b0
    lowerCI_b0<<-b0-t_table*s_b0

    if(B10==0){
      #t-test
      t_star<<-(b0-B10)/s_b0
      t_table<<-qt(1-alpha/2,n-p)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_b0 & B10>=lowerCI_b0) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-2*pt(-abs(t_star),n-p)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }
    else{
      #t-test
      t_star<<-(b0-B10)/s_b0
      t_table<<-qt(1-alpha,n-p)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_b0 & B10>=lowerCI_b0) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-pt(-abs(t_star),n-p)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }

    t_star_b0<<-t_star
    t_table_b0<<-t_table
    p_value_b0<<-p_value


    ################################################

    #inference on B2

    #confidence interval
    if(multivariable==TRUE){
      upperCI_b2<<-b2+t_table*s_b2
      lowerCI_b2<<-b2-t_table*s_b2

      if(B10==0){
        #t-test
        t_star<<-(b2-B10)/s_b2
        t_table<<-qt(1-alpha/2,n-p)

        if(abs(t_star)<=t_table) print("conclude H0")
        else print("conclude Ha")

        if(B10<=upperCI_b2 & B10>=lowerCI_b2) print("conclude H0 because B10 lies inside CI")
        else print("conclude Ha because B10 does not lie inside CI")

        p_value<<-2*pt(-abs(t_star),n-p)
        if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
        else print("conclude Ha because p_value is smaller than significance level aka. alpha")
      }
      else{
        #t-test
        t_star<<-(b2-B10)/s_b2
        t_table<<-qt(1-alpha,n-p)

        if(abs(t_star)<=t_table) print("conclude H0")
        else print("conclude Ha")

        if(B10<=upperCI_b2 & B10>=lowerCI_b2) print("conclude H0 because B10 lies inside CI")
        else print("conclude Ha because B10 does not lie inside CI")

        p_value<<-pt(-abs(t_star),n-p)
        if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
        else print("conclude Ha because p_value is smaller than significance level aka. alpha")
      }
    }


    t_star_b2<<-t_star
    t_table_b2<<-t_table
    p_value_b2<<-p_value


    ################################################

    #Inference on Yh_hat (interval estimation)
    Yh_hat<<-b0+b1*Xh1
    if(multivariable==TRUE)Yh_hat<<-b0+b1*Xh1+b2*Xh2

    #estimated variance
    s2_Yh_hat<<-MSE*t(Xh)%*%solve(t(X)%*%X)%*%Xh
    s_Yh_hat<<-sqrt(s2_Yh_hat)

    #confidence interval
    upperCI_Yh_hat<<-Yh_hat+t_table*s_Yh_hat
    lowerCI_Yh_hat<<-Yh_hat-t_table*s_Yh_hat


    if(B10==0){
      #t-test
      t_star<<-(Yh_hat-B10)/Yh_hat
      t_table<<-qt(1-alpha/2,n-p)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_Yh_hat & B10>=lowerCI_Yh_hat) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-2*pt(-abs(t_star),n-p)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }
    else {
      #t-test
      t_star<<-(Yh_hat-B10)/Yh_hat
      t_table<<-qt(1-alpha,n-2)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_Yh_hat & B10>=lowerCI_Yh_hat) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-pt(-abs(t_star),n-2)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }
    t_star_Yh_hat<<-t_star
    t_table_Yh_hat<<-t_table
    p_value_Yh_hat<<-p_value


    ################################################

    #Inference on pred (prediction of new observation estimation for given Xh)
    #estimated variance
    s2_pred<<-MSE*(1+MSE*t(Xh)%*%solve(t(X)%*%X)%*%Xh)
    s_pred<<-sqrt(s2_pred)


    #prediction interval
    upperPI_pred<<-Yh_hat+t_table*s_pred
    lowerPI_pred<<-Yh_hat-t_table*s_pred

    ################################################

    #Inference on predmean (prediction mean of m new observation for given Xh)
    #estimated variance
    s2_predmean<<-MSE*(1/m+MSE*t(Xh)%*%solve(t(X)%*%X)%*%Xh)
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

    #Upper_emperical_rule_68<<-sqrt(MSE)
    #lower_emperical_rule_68<<-(-sqrt(MSE))

    #Upper_emperical_rule_90<<-1.645*sqrt(MSE)
    #lower_emperical_rule_90<<-(-1.645*sqrt(MSE))

    #count1<<-0
    #for (i in residual){
    #  if(residual[i]<Upper_emperical_rule_68 & residual[i]>lower_emperical_rule_68) count1<<-count1+1
    #}

    #count2<<-0
    #for (i in residual){
    #  if(residual[i]<Upper_emperical_rule_90 & residual[i]>lower_emperical_rule_90) count2<<-count2+1
    #}

    #count1<<-count1/n*100
    #count2<<-count2/n*100

    #if(count1>=68 & count2>=90) print("Emperical rule respected")
    #else print("Emperical rule not respected")

    plot(sort(Ee),sort(residual))

    fit2<- lm(sort(residual) ~ sort(Ee))
    anova(fit2)
    summary(fit2)

    ################################################

    #F_test_for_lack_of_fit



    ################################################

    #bonferroni_inequality_to_estimate_B0_and_B1
    B_CI<<-qt(1-alpha/4,n-p)

    upperCI_bonferroni_B1<<-b1+B_CI*s_b1
    lowerCI_bonferroni_B1<<-b1-B_CI*s_b1

    upperCI_bonferroni_B0<<-b0+B_CI*s_b0
    lowerCI_bonferroni_B0<<-b0-B_CI*s_b0

    ################################################

    #working_hotteling_procedure_for_simultaneous_estimation_of_mean_response
    W2<<-2*qf(1-alpha,2,n-p)
    W<<-sqrt(W2)

    upperCI_WH_Yh_hat<<-Yh_hat+W*s_Yh_hat
    lowerCI_WH_Yh_hat<<-Yh_hat-W*s_Yh_hat

    ################################################

    #bonferroni_procedure_for_simultaneous_estimation_of_mean_response
    B_CP<<-qt(1-(alpha/(2*g)),n-p)

    upperCP_bonferroni_Yh_hat<<-Yh_hat+B_CP*s_Yh_hat
    lowerCP_bonferroni_Yh_hat<<-Yh_hat-B_CP*s_Yh_hat

    ################################################

    #scheffe_procedure_for_simultaneous_PI_for_g_new_observations
    S2<<-2*qf(1-alpha,g,n-p)
    S<<-sqrt(S2)

    upperPI_scheffe_pred<<-Yh_hat+S*s_pred
    lowerPI_scheffe_pred<<-Yh_hat-S*s_pred

    ################################################

    #bonferroni_procedure_for_simultaneous_PI_for_g_new_observations
    upperPI_bonferroni_pred<<-Yh_hat+B_CP*s_pred
    lowerPI_bonferroni_pred<<-Yh_hat-B_CP*s_pred
  }


















  else{
    #variables
    X<<- data1$X1var
    Y<<- data1$Yvar
    n<<-length(X)
    fit <<-lm(Y~X)
    print(summary(lm(Y~X)))

    #variables assignment
    alpha<<-num1
    B10<<-num2
    Xh<<-num3
    m<<-num4
    g<<-num5

    #calculations
    Sxy<<-sum((X-mean(X))*(Y-mean(Y)))
    Sxx<<-sum((X-mean(X))^2)
    Syy<<-sum((Y-mean(Y))^2)

    #regression line coefficients
    b1<<-Sxy/Sxx
    b0<<-mean(Y)-b1*mean(X)
    Fitted<<-b0+b1*X

    #residuals
    residual<<-Y-Fitted

    SSE<<-sum(residual^2)
    MSE<<-SSE/(n-2)

    plot(X,Y)
    lines(X,Fitted)

    plot(sort(X),residual)
    lines(sort(X),rep(0,n))

    #estimation of residual variance sigma_Sq
    s2<<- MSE
    s<<- sqrt(MSE)

    #################################################

    #ANOVA

    #semistudentized residuals
    e_star1<<-(residual-mean(residual))/sqrt(MSE)

    #prediction deviation (error)
    prediction_deviation<<-Fitted-mean(Y)

    #total deviation
    total_deviation<<-Y-mean(Y)

    #sum of squares
    SSTO<<-sum(total_deviation^2)
    SSR<<-sum(prediction_deviation^2)
    SSE<<-sum(residual^2)

    #mean of squares
    MSTO<<-SSTO/(n-1)
    MSR<<-SSR/1
    MSE<<-SSE/(n-2)

    #F-test for linear regression
    F_star<<-MSR/MSE
    F_table<<-qf(1-alpha,1,n-2)
    if(F_star<=F_table) print("conclude H0")
    else print("conclude Ha")
    print(paste0("for linear regression: F_star: ",F_star," F_table: ",F_table))

    #t-test (using F) for linear regression
    t_star_using_F<<-sqrt(F_star)
    t_table_using_F<<-qt(1-alpha/2,n-2)
    if(t_star_using_F<=t_table_using_F) print("conclude H0")
    else print("conclude Ha")
    print(paste0("t_star_using_F: ",t_star_using_F," t_table_using_F: ",t_table_using_F))

    #ANOVA in R
    print(anova(fit))

    #coefficient of determination (R^2) and correlation coefficient (r)
    R_2<<-SSR/SSTO
    r<<-sqrt(R_2)
    if (b1>0) print(paste0("r: +",r))
    else print(paste0("r: -",r))

    ################################################

    #inference on B1

    #estimated variance
    s2_b1<<-MSE/Sxx
    s_b1<<-sqrt(s2_b1)

    #t-test concerning B1
    if(B10==0){
      t_star<<-(b1-B10)/s_b1
      t_table<<-qt(1-alpha/2,n-2)

      p_value<<-2*pt(-abs(t_star),n-2)
    }
    else{
      t_star<<-(b1-B10)/s_b1
      t_table<<-qt(1-alpha,n-2)

      p_value<<-pt(-abs(t_star),n-2)
    }

    #confidence interval of B1
    upperCI_b1<<-b1+t_table*s_b1
    lowerCI_b1<<-b1-t_table*s_b1



    if(abs(t_star)<=t_table) print("conclude H0")
    else print("conclude Ha")

    if(B10<=upperCI_b1 & B10>=lowerCI_b1) print("conclude H0 because B10 lies inside CI")
    else print("conclude Ha because B10 does not lie inside CI")

    if(p_value>alpha) print("conclude H0 because p_value is larger than significance level: alpha")
    else print ("conclude Ha because p_value is smaller than significance level: alpha")

    t_star_b1<<-t_star
    t_table_b1<<-t_table
    p_value_b1<<-p_value
    ################################################

    #inference on B0

    #estimated variance
    s2_b0<<-MSE*(1/n+(mean(X))^2/Sxx)
    s_b0<<-sqrt(s2_b0)

    #confidence interval
    upperCI_b0<<-b0+t_table*s_b0
    lowerCI_b0<<-b0-t_table*s_b0


    if(B10==0){
      #t-test
      t_star<<-(b0-B10)/s_b0
      t_table<<-qt(1-alpha/2,n-2)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_b0 & B10>=lowerCI_b0) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-2*pt(-abs(t_star),n-2)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }
    else{
      #t-test
      t_star<<-(b0-B10)/s_b0
      t_table<<-qt(1-alpha,n-2)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_b0 & B10>=lowerCI_b0) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-pt(-abs(t_star),n-2)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }

    t_star_b0<<-t_star
    t_table_b0<<-t_table
    p_value_b0<<-p_value

    ################################################

    #Inference on Yh_hat (interval estimation)
    Yh_hat<<-b0+b1*Xh

    #estimated variance
    s2_Yh_hat<<-MSE*(1/n+(Xh-mean(X))^2/Sxx)
    s_Yh_hat<<-sqrt(s2_Yh_hat)

    #confidence interval
    upperCI_Yh_hat<<-Yh_hat+t_table*s_Yh_hat
    lowerCI_Yh_hat<<-Yh_hat-t_table*s_Yh_hat

    if(B10==0){
      #t-test
      t_star<<-(Yh_hat-B10)/Yh_hat
      t_table<<-qt(1-alpha/2,n-2)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_Yh_hat & B10>=lowerCI_Yh_hat) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-2*pt(-abs(t_star),n-2)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }
    else {
      #t-test
      t_star<<-(Yh_hat-B10)/Yh_hat
      t_table<<-qt(1-alpha,n-2)

      if(abs(t_star)<=t_table) print("conclude H0")
      else print("conclude Ha")

      if(B10<=upperCI_Yh_hat & B10>=lowerCI_Yh_hat) print("conclude H0 because B10 lies inside CI")
      else print("conclude Ha because B10 does not lie inside CI")

      p_value<<-pt(-abs(t_star),n-2)
      if(p_value>alpha) print("conclude H0 because p_value is larger than significance level aka. alpha")
      else print("conclude Ha because p_value is smaller than significance level aka. alpha")
    }

    t_star_Yh_hat<<-t_star
    t_table_Yh_hat<<-t_table
    p_value_Yh_hat<<-p_value

    ################################################

    #Inference on pred (prediction of new observation estimation for given Xh)
    #estimated variance
    s2_pred<<-MSE*(1+1/n+(Xh-mean(X))^2/Sxx)
    s_pred<<-sqrt(s2_pred)


    #prediction interval
    upperPI_pred<<-Yh_hat+t_table*s_pred
    lowerPI_pred<<-Yh_hat-t_table*s_pred

    ################################################

    #Inference on predmean (prediction mean of m new observation for given Xh)
    #estimated variance
    s2_predmean<<-MSE*(1/m+1/n+(Xh-mean(X))^2/Sxx)
    s_predmean<<-sqrt(s2_predmean)

    #prediction interval
    upperPI_predmean<<-Yh_hat+t_table*s_predmean
    lowerPI_predmean<<-Yh_hat-t_table*s_predmean

    ################################################

    #Prediction band
    s2_Yh_hat<<-MSE*(1/n+(Xh-mean(X))^2/Sxx)
    s_Yh_hat<<-sqrt(s2_Yh_hat)

    W2<<-2*qf(1-alpha,2,n-2)
    W<<-sqrt(W2)

    upperPB<<-Yh_hat+W*s_Yh_hat
    lowerPB<<-Yh_hat-W*s_Yh_hat

    ################################################

    # general test approach
    SSE_F_general_test<<-SSE
    SSE_R_general_test<<-SSTO

    #F test (test statistic of {SSE(R) - SSE(F)} )
    df_F_general_test<<-n-2
    df_R_general_test<<-n-1

    F_star_general_test<<-((SSE_R_general_test-SSE_F_general_test)/
                             (df_R_general_test-df_F_general_test))/
                                ((SSE_F_general_test)/df_F_general_test)

    F_table_general_test<<-qf(1-alpha, df_R_general_test-df_F_general_test, df_F_general_test)

    if(F_star_general_test<=F_table_general_test) print("general test approach: conclude H0")
    else print("conclude Ha")

    ################################################

    #exepected_value_of_e_under_normality
    k=1:n
    Ee<<-sqrt(MSE)*qnorm((k-0.375)/(n+0.25))

    #Upper_emperical_rule_68<<-sqrt(MSE)
    #lower_emperical_rule_68<<-(-sqrt(MSE))

    #Upper_emperical_rule_90<<-1.645*sqrt(MSE)
    #lower_emperical_rule_90<<-(-1.645*sqrt(MSE))

    #count1<<-0
    #for (i in residual){
    #  if(residual[i]<Upper_emperical_rule_68 & residual[i]>lower_emperical_rule_68) count1<<-count1+1
    #}

    #count2<<-0
    #for (i in residual){
    #  if(residual[i]<Upper_emperical_rule_90 & residual[i]>lower_emperical_rule_90) count2<<-count2+1
    #}

    #count1<<-count1/n*100
    #count2<<-count2/n*100

    #if(count1>=68 & count2>=90) print("Emperical rule respected")
    #else print("Emperical rule not respected")

    plot(sort(Ee),sort(residual))
    print("FIT2")
    fit2<- lm(sort(residual) ~ sort(Ee))
    # print(anova(fit2))
    # print(summary(fit2))

    ################################################

    #F_test_for_lack_of_fit
    data_table<<-data.frame(x=X,y=Y,nj=c(rep(1,n)),
                          mean=c(rep(0,n)),stDev=c(rep(NA,n)))

    for(i in 1:nrow(data_table)){
      if(data_table[i,1]==0)next
      for(j in i:nrow(data_table)){
        if(i==j) next
        if(data_table[j,1]==data_table[i,1]) {

          data_table[i,2]<<-paste0(data_table[i,2],", ",data_table[j,2])
          data_table[i,3]<<-data_table[i,3]+1

          data_table[j,]<<-0
        }
      }
    }

    data_table<<-subset(data_table,x!=0)
    rownames(data_table) <<- NULL

    c_<<-nrow(data_table)

    for(i in 1:nrow(data_table)){
      v_=c()
      for(j in 1:data_table[i,3]){
        v_=cbind(v_,as.double(tokenize_string(data_table[i,2])
                              [1+(j-1)*3,1]))
        data_table[i,4]<<-mean(v_)
        data_table[i,5]<<-sd(v_)
      }
    }

    data_table[is.na(data_table)]<<-0
    print(data_table)

    SSPE<<-0
    for(i in 1:nrow(data_table)){
      SSPE<<-SSPE+data_table[i,5]^2*(data_table[i,3]-1)
    }

    MSPE<<-SSPE/(n-c_)
    SSLF<<-(SSE-SSPE)
    MSLF<<-SSLF/(c_-2)

    F_star_LF<<-MSLF/MSPE
    F_table_LF<<-qf(1-alpha, c_-2, n-c_)

    if(F_star<=F_table) print("conclude H0")
    else print("conclude Ha")

    Full<<-lm(Y~0 + as.factor(X))
    print(anova(fit,Full))

    ################################################

    #bonferroni_inequality_to_estimate_B0_and_B1
    B_CI<<-qt(1-alpha/4,n-2)

    upperCI_bonferroni_B1<<-b1+B_CI*s_b1
    lowerCI_bonferroni_B1<<-b1-B_CI*s_b1

    upperCI_bonferroni_B0<<-b0+B_CI*s_b0
    lowerCI_bonferroni_B0<<-b0-B_CI*s_b0

    ################################################

    #working_hotteling_procedure_for_simultaneous_estimation_of_mean_response
    W2<<-2*qf(1-alpha,2,n-2)
    W<<-sqrt(W2)

    s2_Yh_hat<<-MSE*(1/n+(Xh-mean(X))^2/Sxx)
    s_Yh_hat<<-sqrt(s2_Yh_hat)

    upperCI_WH_Yh_hat<<-Yh_hat+W*s_Yh_hat
    lowerCI_WH_Yh_hat<<-Yh_hat-W*s_Yh_hat

    ################################################

    #bonferroni_procedure_for_simultaneous_estimation_of_mean_response
    B_CP<<-qt(1-(alpha/(2*g)),n-2)

    s2_Yh_hat<<-MSE*(1/n+(Xh-mean(X))^2/Sxx)
    s_Yh_hat<<-sqrt(s2_Yh_hat)

    upperCP_bonferroni_Yh_hat<<-Yh_hat+B_CP*s_Yh_hat
    lowerCP_bonferroni_Yh_hat<<-Yh_hat-B_CP*s_Yh_hat

    ################################################

    #scheffe_procedure_for_simultaneous_PI_for_g_new_observations
    S2<<-2*qf(1-alpha,g,n-2)
    S<<-sqrt(S2)

    upperPI_scheffe_pred<<-Yh_hat+S*s_pred
    lowerPI_scheffe_pred<<-Yh_hat-S*s_pred

    ################################################

    #bonferroni_procedure_for_simultaneous_PI_for_g_new_observations
    upperPI_bonferroni_pred<<-Yh_hat+B_CP*s_pred
    lowerPI_bonferroni_pred<<-Yh_hat-B_CP*s_pred
  }

}

# statistic_test <-function(test, point_estimator, sd_pt_estim, num1, num2){
#   if(test==t){
#     if(num2==0){
#       t_star<<-(point_estimator-B10)/sd_pt_estim
#       t_table<<-qt(1-alpha/2,n-2)
#
#       p_value<<-2*pt(-abs(point_estimator),n-2)
#     }
#     else{
#       t_star<<-(point_estimator-B10)/s_b1
#       t_table<<-qt(1-alpha,n-2)
#
#       p_value<<-pt(-abs(test_statistic),n-2)
#     }
#   }
#
#   else if(test==F){
#     F_star<<- MSR/MSE
#     F_table<<- qf(1-alpha,1,n-2)
#   }
#   else if(test==t_using_F){
#     t_star<<-sqrt(F_star)
#     t_table<<-qt(1-alpha/2,n-2)
#   }
#   else if(test==F_general){
#     df_F_general_test<<-n-2
#     df_R_general_test<<-n-1
#
#     F_star_general_test<<-((SSE_R-SSE_F)/(df_R-df_F))/((SSE_F)/df_F)
#
#     F_table_general_test<<-qf(1-alpha, df_R-df_F, df_F)
#   }
#   else if(test==F_LF){
#     F_star_LF<<-MSLF/MSPE
#     F_table_LF<<-qf(1-alpha, c-2, n-c)
#   }
#   #confidence interval
#   upper<<-point_estimator+value_table*sd_pt_estim
#   lower<<-point_estimator-value_table*sd_pt_estim
# }
#
# test_conclusion <- function(B10, statistic_test, value_table, lower, upper){
#   if(abs(statistic_test)<=value_table) print("conclude H0")
#   else print("conclude Ha")
#
#   if(B10<=upper & B10>=lower) print("conclude H0 because B10 lies inside CI")
#   else print("conclude Ha because B10 does not lie inside CI")
#
#   if(p_value>alpha) print("conclude H0 because p_value is larger than significance level: alpha")
#   else print ("conclude Ha because p_value is smaller than significance level: alpha")
# }



