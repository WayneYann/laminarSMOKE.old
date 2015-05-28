#if OPENSMOKE_USE_ISAT == 1
if(isatCheck == true)
{
	#include "chemistry_ISAT.H"
}
else
{
	if (drg_analysis == false)
		#include "chemistry_DI.H"
	else
		#include "chemistry_DRG.H"
}
#else
	if (drg_analysis == false)
		#include "chemistry_DI.H"
	else
		#include "chemistry_DRG.H"
#endif
