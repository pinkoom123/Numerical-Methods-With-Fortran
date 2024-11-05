PROGRAM hubble_constant
   IMPLICIT NONE
   
   REAL, DIMENSION(15) :: d, s  ! Arrays for distances and speeds
   REAL, DIMENSION(2,2) :: A, A_inv  ! 2x2 matrix 
   REAL, DIMENSION(2) :: B, X  
   REAL :: sum_d, sum_s, sum_ds, sum_d2, detA,mean_d
   INTEGER :: ios, i, N
   REAL :: ss_total, ss_residual, R2, se_beta1, se_beta2, mean_s, residual_var
 
   N = 15

   ! Open the file to read the data
   OPEN(30, FILE="hubble_data.txt", STATUS='old', ACTION='read', IOSTAT=ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error opening file'
      STOP
   END IF

   ! Read the distances (d) and speeds (s) from the file
   READ(30,*) ! Skip header
   DO i = 1, N
      READ(30,*) s(i), d(i)
   END DO
       !PRINT*, s
       !PRINT*, d
   
   CLOSE(30)
   
    
   ! Initialize sums for the linear fit
   sum_d = SUM(d)
   sum_s = SUM(s)
   sum_ds = SUM(d * s)  ! Sum of products of distances and speeds
   sum_d2 = SUM(d**2)   ! Sum of squares of distances

   ! Fill matrix A 
   A(1,1) = N
   A(1,2) = sum_d
   A(2,1) = sum_d
   A(2,2) = sum_d2
   
   ! print*, A
   ! Fill vector B 
   B(1) = sum_s
   B(2) = sum_ds

   ! Calculate the determinant of matrix A
   detA = A(1,1) * A(2,2) - A(1,2) * A(2,1)

   ! Calculate the inverse of matrix A
   A_inv(1,1) =  A(2,2) / detA
   A_inv(1,2) = -A(1,2) / detA
   A_inv(2,1) = -A(2,1) / detA
   A_inv(2,2) =  A(1,1) / detA
   
    !PRINT*,A_inv
   
   ! Solve for X (beta1, beta2)
   X(1) = A_inv(1,1) * B(1) + A_inv(1,2) * B(2)
   X(2) = A_inv(2,1) * B(1) + A_inv(2,2) * B(2)
  
    ! ANOTHER WAY TO CALCULATE THE ESTIMATES
    ! X = MATMUL(A_INV,B)
    
     ! Output the results
   !PRINT*, 'The intercept (beta1) is: ', X(1)
   !PRINT*, 'The slope (beta2, Hubble constant) is: ', X(2)
   
   ! Calculate the coefficient of determination (R^2)
   mean_s = sum_s / N
   ss_total = SUM((s - mean_s)**2)  ! Total sum of squares (variance in data)
   ss_residual = SUM((s - (X(1) + X(2) * d))**2)  ! Residual sum of squares
   R2 = 1 - (ss_residual / ss_total)

   print*,ss_residual
   ! Calculate the variance squared
   residual_var = ss_residual / (N - 2)
   PRINT*,'The variance squared  :',residual_var
   
   ! Calculate the mean of distances
   mean_d = sum_d / N

   ! Calculate standard errors for beta1 and beta2 using the correct formulas
   se_beta1 = SQRT(residual_var * (1.0/N + mean_d**2 / SUM((d - mean_d)**2)))
   se_beta2 = SQRT(residual_var / SUM((d - mean_d)**2))
   

! Output the results to a text file
   OPEN(40, FILE="fit.txt", STATUS='replace', ACTION='write', IOSTAT=ios)
   IF (ios /= 0) THEN
      PRINT*, 'Error opening output file'
      STOP
   END IF

   WRITE(40,*) 'Estimates for the linear fit v_r(d) = beta1 + beta2 * d'
   WRITE(40,*) '---------------------------------------------'
   WRITE(40,*) 'Intercept (beta1)          : ', X(1)
   WRITE(40,*) 'Slope (beta2, Hubble constant) : ', X(2)
   WRITE(40,*) '---------------------------------------------'
   WRITE(40,*) 'Standard error for beta1    : ', se_beta1
   WRITE(40,*) 'Standard error for beta2    : ', se_beta2
   WRITE(40,*) 'Coefficient of determination (R^2) : ', R2

   CLOSE(40)

   ! Output results to console as well
   PRINT*, 'The intercept (beta1) is: ', X(1)
   PRINT*, 'The slope (beta2, Hubble constant) is: ', X(2)
   PRINT*, 'R^2 (coefficient of determination) is: ', R2


   
END PROGRAM hubble_constant



