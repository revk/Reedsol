// $Log: reedsol.h,v $
// Revision 1.3  2006/01/25 11:10:43  cvs
// New sqlhtml
//
// This software is provided under the terms of the GPL v2 or later.
// This software is provided free of charge with a full "Money back" guarantee.
// Use entirely at your own risk. We accept no liability. If you don't like that - don't use it.

// Revision 1.2  2004/09/09 07:45:09  cvs
// Added change history to source files
// Added "info" type to IEC16022
// Added exact size checking shortcodes on encoding generation for iec16022
//

//#define	RS_CONNECT	// If correction is needed too

void rs_init(int poly,int nsym, int index);
void rs_encode(int len, unsigned char *data,unsigned char*res);
#ifdef RS_CORRECT
int rs_correct(int datalen, byte *data);
#endif

