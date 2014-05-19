16 May 2004 immoments added
            ploto     added
            featsel routines with crossval implemented
25 Jun 2004 fdsc      added
30 Aug 2004 proxm     homogeneous measure added (x*y)^p
            disperror, confmat: writing to files enabled
26 Sep 2004 gaussm, emclust: get right number of components
            crisp / soft label handling for these routines
30 Sep 2004 plotr     allow multiple error curves
            scalem    variance option corrected
21 Oct 2004 svc_nu, svo_nu, added
 6 Dec 2004 parzendc, parzen_map: soft labels
            nmc: bug fix, allow soft labels
            prtest: bug fix, label generation
12 Dec 2004 plotr: cell array of error structures
14 Dec 2004 setlabels, roc, normm: bug fixes to correct / prevent
            reordering of lablist
20 Jan 2005 prwarning level now default 1
            some warning levels adapted
            in particular for empty priors
 1 Feb 2005 clevalf, debugged from reusing the same testset
24 Apr 2005 labcmp added
29 Oct 2005 trying to avoid prior warnings
  
----------------------------

 8 Mar 2006 start datafile options
26 Sep 2006 multiple labeling, see multi-labeling
26 Sep 2006 structure field handling for user and ident field in datasets
28 Nov 2006 testauc added

18 Dec 2006 global variables like PRMEMORY, PRTRACE, GRIDSIZE, etc 
            are made persistent instead of global
19 Dec 2006 plotr renamed in plote (make space for a regression plot command)
            svc svc_nu: general option for kernel computation
22 Mar 2007 datafile construct ready
19 Apr 2007 Normal_map: inversion at learning
            Gaussm, mogs, improved: noise, more trials, regularisation
17 May 2007 NAIVEBC generalised with PARZENM and GAUSSM
18 May 2007 Streamlining PRPROGRESS and CLOSEMESS
            NMC now based on normal distributions with identically shaped
            spherical classifiers. This gives a better soft output. Consequently
            the soft version of EMCLUST is improved as it uses NMC for 
            initialisation.
 5 Jun 2007 TESTC extended with various criteria
            Automatic parameter optimization in some routines
            RSSCC and RBSVC added
 3 Jul 2008 PRPROGRESS partially replaced by PRWAITBAR
10 Jul 2008 ident reorganised
21 Jul 2008 createdatafile
18 Aug 2008 Image Database Retrieval GUI
 5 Jan 2009 im_patch, band2obj and bandsel added
- new trainable combiners: naivebcc and mlrc (thanks to Chunxia!)
- Leave One Set Out crossvalidation: loso
- testc accepts missing classes in case of no_priors
- implementation of multi-band time-signals as multiband 1-D images
- this script generates now proper p-code :-)
10 Mar 2009 DCSC added
18 Mar 2009 MODSELC added

April  2011 FDSC, MDSC, DISNORM added
