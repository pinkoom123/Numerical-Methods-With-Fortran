MODULE N_linear
   IMPLICIT NONE
          CONTAINS
        
    ! Function for Planck's law model B(lamda; beta1, beta2)
    REAL FUNCTION beta(w, b1, b2)
          IMPLICIT NONE
          REAL, INTENT(IN) :: w, b1, b2
          REAL :: exp_term
          exp_term = EXP(b2 / w)
          beta = b1 / (w**5 * (exp_term - 1.0))
       END FUNCTION beta

    ! Partial derivative with respect to beta1
    REAL FUNCTION Derivative_b1(w, b2)
           IMPLICIT NONE
           REAL, INTENT(IN) :: w, b2
           REAL :: exp_term
           exp_term = EXP(b2 / w)
           Derivative_b1 = 1.0 / (w**5 * (exp_term - 1.0))
    END FUNCTION Derivative_b1

    ! Partial derivative with respect to beta2
   REAL FUNCTION Derivative_b2(w, b1, b2)
            IMPLICIT NONE
            REAL, INTENT(IN) :: w, b1, b2
            REAL :: exp_term
            exp_term = EXP(b2 / w)
            Derivative_b2 = -b1 * exp_term / (w**6 * (exp_term - 1.0)**2)
    END FUNCTION Derivative_b2
      
  
END MODULE N_linear


PROGRAM Planks_law
  USE N_linear
 IMPLICIT NONE  
   REAL,DIMENSION(1612):: w,b,residuals
   REAL:: alpha,e,sum_res,DET,var_error,T,constant
   REAL,DIMENSION(2):: beta_vals,new_beta,delta_beta,SE_beta
   REAL,DIMENSION(1612,2)::J
   REAL,DIMENSION(2,2):: JTJ_inv,JTJ,var_beta
   INTEGER:: ios,i,it,maxIter,N
   
   N=1612
   constant =14387.77
! initialize values
    alpha = 0.1
    maxIter=1000
    e= 1.0E-6
    beta_vals = (/1.0,1.0/) !intial guess

! Open the file to read the data
   OPEN(30, FILE="sun_data.txt", STATUS='old', ACTION='read', IOSTAT=ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error opening file'
      STOP
   END IF

   ! Read the wavelength (w) and intensity (b) from the file
   READ(30,*) ! Skip header
    READ(30,*) ! Skip header
   DO i = 1, 1612
      READ(30,*) w(i), b(i)
   END DO
       !PRINT*, w
       !PRINT*, b
   
   CLOSE(30)
   
 ! iterative Gauss-Newton method
    DO it= 1, maxIter
        ! Compute the residuals and Jacobian matrix J
        sum_res = 0.0
        DO i = 1, 1612
            residuals(i) = b(i) - beta(w(i), beta_vals(1), beta_vals(2))
            sum_res = sum_res + residuals(i)**2

            J(i, 1) = Derivative_b1(w(i), beta_vals(2))
            J(i, 2) = Derivative_b2(w(i), beta_vals(1), beta_vals(2))
 
        END DO

     ! Compute JTJ and JTr (transpose of Jacobian times residuals)
     JTJ = MATMUL(TRANSPOSE(J), J)  ! JTJ = J^T * J
     new_beta = MATMUL(TRANSPOSE(J),residuals)

     ! Compute the determinant of JTJ
        det = JTJ(1,1) * JTJ(2,2) - JTJ(1,2) * JTJ(2,1)

       
        ! Compute the inverse of JTJ
        JTJ_inv(1,1) =  JTJ(2,2) / det
        JTJ_inv(1,2) = -JTJ(1,2) / det
        JTJ_inv(2,1) = -JTJ(2,1) / det
        JTJ_inv(2,2) =  JTJ(1,1) / det

        ! Solve for delta_beta: delta_beta = inv(J^T * J) * J^T * r
        delta_beta = MATMUL(JTJ_inv, new_beta)
        
        ! Check for convergence:
        IF (SUM(delta_beta**2) <= e) EXIT
    
        ! Update beta values with step size alpha
        beta_vals = beta_vals + alpha * delta_beta
    END DO
     ! Calculate the variance of the error (residual variance)
    var_error = sum_res / (N - 2)

    ! Calculate the covariance matrix for beta (var_beta)
    var_beta = var_error * JTJ_inv  ! 

    ! Compute the standard errors: SE(beta_i) = sqrt(Var(beta_i))
    SE_beta(1) = SQRT(var_beta(1,1))
    SE_beta(2) = SQRT(var_beta(2,2))
    
    ! calculate for temperature
      T = constant/beta_vals(2)
 
  
  ! Output the results to a text file
   OPEN(40, FILE="fitt.txt", STATUS='replace', ACTION='write', IOSTAT=ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error opening output file'
      STOP
   END IF

   WRITE(40,*) 'Estimates for the non-linear fit β ← β + α∆β.'
   WRITE(40,*) '---------------------------------------------'
   WRITE(40,*) 'Estimate (beta1) : ', beta_vals(1)
   WRITE(40,*) 'Estimate (beta2) : ', beta_vals(2)
   WRITE(40,*) '---------------------------------------------'
   WRITE(40,*) 'Standard error for beta1    : ', SE_beta(1)
   WRITE(40,*) 'Standard error for beta2    : ', SE_beta(2)
   WRITE(40,*) '---------------------------------------------'
   WRITE(40,*) 'The variance of the error  :',var_error
   WRITE(40,*) 'The covariance matrix for beta :',  var_beta
   WRITE(40,*) "The temperature calculated is :" ,T,"Kelvin"
   CLOSE(40)

   ! Output results to console as well
   PRINT*,'The result is written a file fitt.txt'
   PRINT*,""
   PRINT*, 'Estimate (beta1) : ', beta_vals(1)
   PRINT*, 'Estimate (beta2) : ', beta_vals(2)
   PRINT*,""
   PRINT*,'The variance of the error  :',var_error
   PRINT*,"The temperature calculated is :" ,T,"Kelvin"
   PRINT*,""
   PRINT*,'A plot of measured Spectrum vs Fitted Planck Model is stored as a file Non-Linear_graph.png'
   
END PROGRAM Planks_law
