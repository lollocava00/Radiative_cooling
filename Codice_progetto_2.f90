MODULE VARIABILI
    IMPLICIT NONE
	REAL*8:: eccitazione, ionizzazione, ricombinazione, ricombinazione_d, free_free,LAMBDA, &
             n_h0, n_h1, n_he0, n_he1, n_he2, n_e, n_he, n_h,y, T 
    CONTAINS
!TASSI DI IONIZZAZIONE
    REAL*8 FUNCTION gamma_H0()
        gamma_H0=5.85d-11*sqrt(T)/(exp(157809.1d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION gamma_H0
    
    REAL*8 FUNCTION gamma_He0()
        gamma_He0=2.38d-11*sqrt(T)/(exp(285335.4d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION gamma_He0
    
    REAL*8 FUNCTION gamma_He1()
        gamma_He1=5.68d-12*sqrt(T)/(exp(631515.0d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION gamma_He1
!TASSI DI RICOMBINAZIONE
    REAL*8 FUNCTION alpha_H1()
        alpha_H1=8.4d-11/(sqrt(T)*(T/1.0d3)**2.0d-1*(1.0d0+(T/1.0d6)**7.0d-1))
    END FUNCTION alpha_H1

    REAL*8 FUNCTION alpha_d()
        alpha_d=1.9d-3*(1.0d0+0.3d0*exp(-94.0d3/T))/((T**1.5d0)*exp(47.0d4/T))
    END FUNCTION alpha_d
    
    REAL*8 FUNCTION alpha_He1()
        alpha_He1=1.5d-10/(T**0.6353d0)
    END FUNCTION alpha_He1
    
    REAL*8 FUNCTION alpha_He2()
        alpha_He2=3.36d-10/(sqrt(T)*((T/1.0d3)**0.2d0)*(1.0d0+(T/1.0d6)**0.7d0))
    END FUNCTION alpha_He2
!CONTRIBUTI DI RAFFREDDAMENTO [erg/s*cm3] 
    REAL*8 FUNCTION lambda_ex_H0(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_ex_H0=7.50d-19*n*m/(exp(118348.0d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION lambda_ex_H0
    
    REAL*8 FUNCTION lambda_ex_H1(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_ex_H1=5.54d-17*n*m/((T**0.397)*exp(473638.0d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION lambda_ex_H1
    
    REAL*8 FUNCTION lambda_ion_H0(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_ion_H0=1.27d-21*sqrt(T)*n*m/(exp(157809.1d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION lambda_ion_H0
    
    REAL*8 FUNCTION lambda_ion_He0(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_ion_He0=9.38d-22*sqrt(T)*n*m/(exp(285335.4d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION lambda_ion_He0
    
    REAL*8 FUNCTION lambda_ion_He1(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_ion_He1=4.95d-22*sqrt(T)*n*m/(exp(631515.0d0/T)*(1.0d0+sqrt(T/1.0d5)))
    END FUNCTION lambda_ion_He1
    
    REAL*8 FUNCTION lambda_rec_H1(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_rec_H1=8.70d-27*sqrt(T)*n*m/(((T/1.0d3)**0.2d0)*(1.0d0+(T/1.0d6)**0.7d0))
    END FUNCTION lambda_rec_H1
    
    REAL*8 FUNCTION lambda_rec_He1(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_rec_He1=1.55d-26*(T**0.3647d0)*n*m
    END FUNCTION lambda_rec_He1
    
    REAL*8 FUNCTION lambda_rec_He2(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_rec_He2=3.48d-26*sqrt(T)*n*m/(((T/1.0d3)**0.2d0)*(1.0d0+(T/1.0d6)**0.7d0))
    END FUNCTION lambda_rec_He2
    
    REAL*8 function lambda_drec_He1(n,m)
        REAL*8,INTENT(IN):: n, m
        lambda_drec_He1=1.24d-13*(1.0d0+0.3d0/exp(94.0d3/T))*n*m/((T**1.5d0)*exp(47.0d4/T))
    END FUNCTION lambda_drec_He1

    REAL*8 FUNCTION gaunt_ff()
        gaunt_ff=1.1d0+0.34d0/exp(((5.5d0-log(T))**2.0)/3.0d0)
    END FUNCTION gaunt_ff
!cooling function
    REAL*8 FUNCTION lambda_ff(n,n1,n2,n3)
        REAL*8,INTENT(IN):: n1, n2, n3, n
        lambda_ff=1.42d-27*gaunt_ff()*sqrt(T)*(n1+n2+4.0d0*n3)*n
    END FUNCTION lambda_ff

    SUBROUTINE cooling_function()
        REAL*8, DIMENSION(2,2):: matrice1
        REAL*8, DIMENSION(2):: cost1, x1
        REAL*8,DIMENSION(3,3):: matrice2
        REAL*8, DIMENSION(3):: cost2, x2
    !2x2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        matrice1=1.0d0
        cost1=0.0d0
        matrice1(1,2)=-alpha_H1()/gamma_H0()
        cost1(2)=1.0d0
        CALL gauss(matrice1,cost1,2,x1)
        n_h0=x1(1)
        n_h1=x1(2)
    !3x3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        matrice2=1.0d0
        cost2=0.0d0
        matrice2(1,2)=-(alpha_He1()+alpha_d())/gamma_He0()
        matrice2(1,3)=0.0d0
        matrice2(2,1)=0.0d0
        matrice2(2,3)=-alpha_He2()/gamma_He1()
        cost2(3)=y
        CALL gauss(matrice2,cost2,3,x2)
        n_he0=x2(1)
        n_he1=x2(2)
        n_he2=x2(3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n_h0=n_h0*n_h
        n_h1=n_h1*n_h
        n_he0=n_he0*n_h
        n_he1=n_he1*n_h
        n_he2=n_he2*n_h
        n_e=n_h1+n_he1+2.0d0*n_he2
    !contributi al cooling!!!!!!!!!!!!!!!!!!!!!
        eccitazione=lambda_ex_H0(n_e,n_h0)+lambda_ex_H1(n_e,n_he1)
        ionizzazione=lambda_ion_H0(n_e,n_h0)+lambda_ion_He0(n_e,n_he0)+lambda_ion_He1(n_e,n_he1)
        ricombinazione=lambda_rec_H1(n_e,n_h1)+lambda_rec_He1(n_e,n_he1)+lambda_rec_He2(n_e,n_he2)
        ricombinazione_d=lambda_drec_He1(n_e,n_he1)
        free_free=lambda_ff(n_e,n_h1,n_he1,n_he2)
        LAMBDA=eccitazione+ionizzazione+ricombinazione+ricombinazione_d+free_free
        
    END SUBROUTINE cooling_function
    
END MODULE VARIABILI

SUBROUTINE gauss(matrice,cost,n,x)
	INTEGER, INTENT(IN):: n
	REAL*8, INTENT(IN):: matrice(n,n), cost(n)
	REAL*8, INTENT(OUT):: x(n)
	REAL*8:: a(n,n), c(n)
	REAL*8:: fakt, summa
	INTEGER:: i,j,k,l

	a=matrice
	c=cost
	DO i=1,n-1
		DO j=i+1,n
			fakt=a(j,i)/a(i,i)
			DO k=1,n
				a(j,k)=a(j,k)-a(i,k)*fakt
			END DO
			c(j)=c(j)-c(i)*fakt
		END DO
	END DO
	x(n)=c(n)/a(n,n)
	DO i=n-1,1,-1
		summa=c(i)
		DO j=i+1,n
			summa=summa-a(i,j)*x(j)
		END DO
		x(i)=summa/a(i,i)
	END DO
END SUBROUTINE gauss

MODULE INTEGRAZIONE
    USE VARIABILI
    IMPLICIT NONE
    REAL*8, PARAMETER:: passo=1.0d-2 ! costante per calcolare il passo
    REAL*8, PARAMETER:: gamma=5.0d0/3.0d0, k=1.380649d-16, mp=1.67262d-24 ! k e mp sono espresse in cgs
	REAL*8:: u_0,mu_0,lambda_0,T0,time_0,rho_0
    CONTAINS
    !DALLA N.B_2 
    REAL*8 FUNCTION fun(Temp,u)
	    REAL*8, INTENT(in):: Temp, u
	    REAL*8:: mu

        T=Temp
        CALL cooling_function()
        mu=(1.0d0+4.0d0*y)/(1.0d0+y+n_e/n_h)
        fun=Temp-((gamma-1)*u_0*u*mu*mp/k)  
    END FUNCTION fun


    REAL*8 FUNCTION deriv(Temp,u)
	    REAL*8, INTENT(in):: Temp, u
	    REAL*8, PARAMETER:: incr=1.0d-5
        !definizione della derivata della funzione fun dato un h=10**-5
	    deriv=(fun(Temp+incr,u)-fun(Temp,u))/incr
    END FUNCTION deriv

    SUBROUTINE secante(Temp,u)
        IMPLICIT NONE
        REAL*8, PARAMETER:: toll=1.0d-6
	    REAL*8, INTENT(IN):: u
	    REAL*8, INTENT(INOUT):: Temp
	    INTEGER:: conta
	    REAL*8:: delta, T_new

    	conta=0
    	DO
	    	conta=conta+1
            !formula della secante modificata!!!!!!!
            T_new=Temp-fun(Temp,u)/deriv(Temp,u)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    delta=abs(T_new-Temp)
		    IF (delta<toll) THEN
			    EXIT
            ELSE IF (T_new<1.0d4) THEN
			    Temp=1.0d4
			    EXIT
            ELSE IF (conta>100) THEN
			    PRINT*, "la convergenza non è stata raggiunta (secante)."
			    STOP
		    END IF
            Temp=T_new
	    END DO
    END SUBROUTINE secante

    SUBROUTINE dydx(u,f,Temp)
        REAL*8, INTENT(IN):: u
        REAL*8, INTENT(INOUT):: Temp
        REAL*8, INTENT(OUT)::f
        REAL*8:: rho, mu
        
        T=Temp
        CALL cooling_function()
        mu=(1.0d0+4.0d0*y)/(1.0d0+y+n_e/n_h)
        rho=mu*(n_h+y*n_h+n_e)*mp
        f=-rho_0*LAMBDA/(rho*lambda_0)

    END SUBROUTINE dydx
    
    SUBROUTINE Runge_Kutta_2(h,u_old,u_new,Temp)  
        REAL*8, INTENT(IN)::h, u_old
        REAL*8, INTENT(INOUT):: Temp
        REAL*8, INTENT(OUT):: u_new
        REAL*8, PARAMETER:: b_2=2.0d0/3.0d0 ! <-- metodo di Ralston per minimizzare l'errore di troncamento
        REAL*8:: k1,k2,a_11,c_1,b_1
        a_11=1.0d0/(2.0d0*b_2)
        c_1=a_11
        b_1=1.0d0-b_2

        CALL dydx(u_old,k1,Temp)
        u_new=u_old+(3.0d0/4.0d0)*h*k1 
        Temp=Temp+(3.0d0/4.0d0)*h

        CALL dydx(u_new,k2,Temp)
        u_new=u_old+h*(b_1*k1+b_2*k2) !<-- formula finale della secante

        !per ottenere la temperatura al tempo t+dt
        CALL secante(Temp,u_new)
    
    END SUBROUTINE Runge_Kutta_2

END MODULE INTEGRAZIONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM PROGETTO_2
    USE VARIABILI
    USE INTEGRAZIONE
    IMPLICIT NONE
    REAL*8::X,T_step,h
    REAL*8:: u,u_tilde,rho,mu,t_cool,time,Temp
    CHARACTER(30):: file_name
    CHARACTER(5):: str
    INTEGER::i,j
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !STATO DI IONIZZAZIONE E FUNZIONE DI COOLING PER X00.76 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	PRINT*, "inizio il calcolo delle diverse densità e della funzione di cooling per X=0.76"
	PRINT*, " "
	OPEN(1,file="densità_relative.dat")
	OPEN(2,file="cooling_function_X0.76.dat")
	WRITE(1,'(a13,5a18)') "Log10(T[°K])","n_H0/n_htot","n_H+/n_htot", &
		"n_He0/n_Hetot","n_He+/n_Hetot","n_He++/n_Hetot"
	WRITE(2,'(a12,a15,a21,2a20,a18,a25)') "Log10(T[°K])","eccitazione", &
		"ionizzazione","ricombinazione","ric_dielettrica","free_free","totale [erg/s*cm**3]"
    X=0.76d0
    y=(1-X)/(4*X)
    n_h=1.d0
    n_he=y
    T_step=(log10(1.0d8)-log10(1.0d4))/200 ! rappresenta il passo da mantenere
    DO i=0,200
        T=10.0d0**(4.0d0+T_step*i)
        CALL cooling_function()
        WRITE(1,'(f8.5,4x,5f18.11)') log10(T),n_h0,n_h1, &
        n_he0/n_he,n_he1/n_he,n_he2/n_he
        WRITE(2,'(f8.5,6f20.10)') log10(T),log10(eccitazione),log10(ionizzazione), &
        log10(ricombinazione),log10(ricombinazione_d),log10(free_free),log10(lambda)
    END DO
    CLOSE(1)
    CLOSE(2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !FUNZIONE DI COOLING per X=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	PRINT*, "calcolo la funzione di cooling per X=1"
	PRINT*, ""
	OPEN(3,file="cooling_function_X1.0.dat")
	WRITE(3,'(a12,a15,a21,2a20,a18,a25)') "Log10(T[°K])","eccitazione", &
		"ionizzazione","ricombinazione","ric_dielettrica","free_free","totale [erg/s*cm**3]"
    X=1.0d0
    y=(1-X)/(4*X)
	DO i=0,200
		T=10.0d0**(4.0d0+T_step*i)
		CALL cooling_function()
		WRITE(3,'(f8.5,6f20.10)') log10(T),log10(eccitazione),log10(ionizzazione), &
			log10(ricombinazione),ricombinazione_d,log10(free_free),log10(lambda)
    END DO
	CLOSE(3)

	PRINT*, 'i risultati del primo e secondo step sono stati salvati nei file'
    PRINT*, 'densità_relative.dat'
    PRINT*, 'cooling_function_X0.76.dat'
    PRINT*, 'cooling_function_X1.0.dat'
	PRINT*, ''
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !INTEGRAZIONE DELLA ODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PRINT*, "procedo con la risoluzione della ODE per l'evoluzione temporale"
    PRINT*, ''
    X=0.76d0
    y=(1-X)/(4*X)
    !IL CICLO INTEGRERÀ IL SISTEMA PER 3 DIVERSE DENSITÀ DI IDROGENO
    DO i=1,3
        IF(i==1) THEN 
            n_h=0.1d0
            str='_n0.1'
        ELSE IF(i==2) THEN
            n_h=1.d0
            str='_n1.0'
        ELSE
            n_h=10.0d0
            str='_n_10'
        END IF
        n_he=n_h*y
        T0=1.d6
        Temp=T0
        T=T0
        CALL cooling_function()
    !Adimensionalizzazione delle variabili!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !queste sono le condizioni iniziali del sistema
        lambda_0=LAMBDA
        mu_0=(1.0d0+4.0d0*y)/(1.0d0+y+n_e/n_h)
        rho_0=mu_0*(n_h+n_he+n_e)*mp
        u_0=k*Temp/((gamma-1.0d0)*mu_0*mp)
        time_0=u_0*rho_0/lambda_0
    
        u_tilde=1.0d0       !<-- u/u_0
        time=0.0d0
        h=passo !<-- il time-step iniziale che verrà modificato man mano nell'integrazione
        t_cool=0.0d0 !<-- che cambierà durante l'integrazione

        file_name='Temp_evolution'//str //'.dat'
        OPEN(10,file=file_name)
        WRITE(10,"(a18,a20)") "Log10(t/t_cooling)", "Log10(T[°K])"
        !integriamo fino a 10 volte il tempo di cooling iniziale(time_0) 
        DO WHILE(time<time_0*10.0d0)

            CALL Runge_Kutta_2(h,u_tilde,u,Temp)

			WRITE(10,"(f14.7,f20.7)") log10(time/time_0), log10(Temp)

            mu=(1.0d0+4.0d0*y)/(1.0d0+y+n_e/n_h) !<--peso molecolare e densità si evolvono poiche dipendenti da T
			rho=mu*(n_h+n_he+n_e)*mp
            !N.B_3 
			t_cool=u*u_0*rho/LAMBDA
			h=passo*t_cool/time_0
			time=time+h*time_0

            !N.B_1
			IF (Temp==1.d4) THEN
				u_tilde=((1.d0/(gamma-1.d0))*((k*Temp)/(mu*mp)))/u_0
            ELSE
				u_tilde=u
			END IF
		END DO
		CLOSE(10)  
    END DO

    PRINT*, "integrazione terminata, i risultati ottenuti sono stati salvati nei file:"
    PRINT*, 'Temp_evolution_n0.1.dat'
    PRINT*, 'Temp_evolution_n1.0.dat'
    PRINT*, 'Temp_evolution_n_10.dat'

END PROGRAM PROGETTO_2