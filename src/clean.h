
#define	CLEAN_ZONE      10.0  /* zone (in sigma) to consider for processing */
#define CLEAN_STACKSIZE 3000  /* replaces prefs.clean_stacksize  */
                              /* (MEMORY_OBJSTACK in sextractor inputs) */
#define CLEAN_MARGIN    0  /* replaces prefs.cleanmargin which was set based */
                           /* on stuff like apertures and vignet size */

extern int addcleanobj(objstruct *, objliststruct *);
extern int subcleanobj(int, objliststruct *);
extern int clean(int, objliststruct *, objliststruct *, LONG *, double);
