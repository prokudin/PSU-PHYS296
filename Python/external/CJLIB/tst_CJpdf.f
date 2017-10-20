      Program tst_CJpdf
      implicit real*8 (a-h,o-z)
      dimension q(4), pdf(-5:5)
      data q/1.3, 10.,31.6228,100./
      print*,'ISET = '
      read*,iset
      call setCJ(iset)
      open(unit=1,file='tst_CJpdf.out', status='unknown')
      do j=1,4
         write(1,2) q(j)**2
         write(1,5)
 5       format(' x times the labelled pdf')
 2       format(' Q**2=',f10.2,' GeV**2')
         write(1,4)
 4       format(5x,' x',9x,'u',11x,'d',11x,'g',10x,'ub',10x,'db',
     2        10x,'sb',10x,' s',10x,'cb',10x,' c',10x,'bb',10x,' b')
         do l=1,19
            x=l*.05
               do k=-5,5
                  pdf(k)=x*CJpdf(k,x,q(j))
               enddo
            write(1,3) x, pdf(1), pdf(2), pdf(0), pdf(-1), 
     2      pdf(-2),pdf(-3), pdf(3), pdf(-4), pdf(4), pdf(-5), pdf(5)
 3          format(f10.3, 11e12.3)
         enddo
      enddo
      print*, 'Output PDFs in tst_CJpdf.out'
      call exit
      end
