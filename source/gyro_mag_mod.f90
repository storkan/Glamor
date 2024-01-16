MODULE gyro_mag_mod
IMPLICIT NONE
PRIVATE
PUBLIC::magnetic_moments, gyro_array
TYPE magnetic_moments
CHARACTER (len=6):: isotope
DOUBLE PRECISION:: gyrovalue
END TYPE magnetic_moments
!source:
!Journal of Magnetic Resonance 156, 323–326 (2002)
!doi:10.1006/jmre.2002.2554
TYPE (magnetic_moments), DIMENSION(34):: gyro_array
DATA gyro_array(1)%isotope, gyro_array(1)%gyrovalue / '1H',26.7522128d7 /
DATA gyro_array(2)%isotope, gyro_array(2)%gyrovalue / '3H',28.5349779d7/
DATA gyro_array(3)%isotope, gyro_array(3)%gyrovalue / '6Li',3.3971709d7/
DATA gyro_array(4)%isotope, gyro_array(4)%gyrovalue / '7Li',10.3977013d7/
DATA gyro_array(5)%isotope, gyro_array(5)%gyrovalue / '11B',8.5847044d7/
DATA gyro_array(6)%isotope, gyro_array(6)%gyrovalue / '13C',6.728284d7/
DATA gyro_array(7)%isotope, gyro_array(7)%gyrovalue / '3He',-20.3801587d7/
DATA gyro_array(8)%isotope, gyro_array(8)%gyrovalue / '15N',-2.71261804d7/
DATA gyro_array(9)%isotope, gyro_array(9)%gyrovalue / '19F',25.18148d7 /
DATA gyro_array(10)%isotope, gyro_array(10)%gyrovalue / '23Na',7.0808493d7 /
DATA gyro_array(11)%isotope, gyro_array(11)%gyrovalue / '27Al',6.9762715d7 /
DATA gyro_array(12)%isotope, gyro_array(12)%gyrovalue / '29Si',-5.3190d7 /
DATA gyro_array(13)%isotope, gyro_array(13)%gyrovalue / '31P',10.8394d7/
DATA gyro_array(14)%isotope, gyro_array(14)%gyrovalue / '35Cl',2.624198d7/
DATA gyro_array(15)%isotope, gyro_array(15)%gyrovalue / '37Cl',2.184368d7/
DATA gyro_array(16)%isotope, gyro_array(16)%gyrovalue / '45Sc',6.5087973d7/
DATA gyro_array(17)%isotope, gyro_array(17)%gyrovalue / '51V',7.0455117d7/
DATA gyro_array(18)%isotope, gyro_array(18)%gyrovalue / '57Fe',0.8680624d7/
DATA gyro_array(19)%isotope, gyro_array(19)%gyrovalue / '71Ga',8.181171d7/
DATA gyro_array(20)%isotope, gyro_array(20)%gyrovalue / '75As',4.596163d7/
DATA gyro_array(21)%isotope, gyro_array(21)%gyrovalue / '77Se',5.1253857d7/
DATA gyro_array(22)%isotope, gyro_array(22)%gyrovalue / '89Y',-1.3162791d7/
DATA gyro_array(23)%isotope, gyro_array(23)%gyrovalue / '103Rh',0.8468d7/
DATA gyro_array(24)%isotope, gyro_array(24)%gyrovalue / '107Ag',-1.0889181d7/
DATA gyro_array(25)%isotope, gyro_array(25)%gyrovalue / '109Ag',-1.2518634d7/
DATA gyro_array(26)%isotope, gyro_array(26)%gyrovalue / '111Cd',-5.6983131d7/
DATA gyro_array(27)%isotope, gyro_array(27)%gyrovalue / '113Cd',-5.9609155d7/
DATA gyro_array(28)%isotope, gyro_array(28)%gyrovalue / '115Sn',-8.8013d7/
DATA gyro_array(29)%isotope, gyro_array(29)%gyrovalue / '117Sn',-9.58879d7/
DATA gyro_array(30)%isotope, gyro_array(30)%gyrovalue / '119Sn',-10.0317d7/
DATA gyro_array(31)%isotope, gyro_array(31)%gyrovalue / '125Te',-8.5108404d7/
DATA gyro_array(32)%isotope, gyro_array(32)%gyrovalue / '129Xe',-7.452103d7/
DATA gyro_array(33)%isotope, gyro_array(33)%gyrovalue / '183W',1.1282403d7/
DATA gyro_array(34)%isotope, gyro_array(34)%gyrovalue / '187Os',0.6192895d7/
CONTAINS
END MODULE gyro_mag_mod


!PROGRAM shit
! USE gyro_mag_mod !, ONLY: gyro_array, magnetic_moments

!IMPLICIT NONE
!integer::i
!DO i=1, 23
!write(*,*)i, gyro_array(i)%isotope, gyro_array(i)%gyrovalue 
!END do
!END PROGRAM shit
