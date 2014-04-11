
#define	PLIST(ptr, elem)	(((pbliststruct *)(ptr))->elem)

#define	PLISTEXIST(elem)	(plistexist_##elem)

#define	PLISTPIX(ptr, elem)	(*((PIXTYPE *)((ptr)+plistoff_##elem)))

#define	PLISTFLAG(ptr, elem)	(*((FLAGTYPE *)((ptr)+plistoff_##elem)))


/* Dummy type for plist */
typedef	char pliststruct;

typedef struct
{
  int     nextpix;
  int     x, y;
  PIXTYPE value;
} pbliststruct;

/* globals */
int plistexist_value, plistexist_dvalue, plistexist_cdvalue,
  plistexist_flag, plistexist_wflag, plistexist_dthresh, plistexist_var,
  plistoff_value, plistoff_dvalue, plistoff_cdvalue,
  plistoff_flag[MAXFLAG], plistoff_wflag, plistoff_dthresh, plistoff_var,
  plistsize;

void init_plist(PIXTYPE *filter, PIXTYPE *cdwfield);
