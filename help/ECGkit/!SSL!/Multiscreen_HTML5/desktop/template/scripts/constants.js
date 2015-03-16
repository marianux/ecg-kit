var PAGENODE = "item";
var BOOKNODE = "book";
var URLNODE = "remoteitem";
var PROJNODE = "project";
var INDEXNODE = "index";
var DATANODE = "data";
var CHUNKINFONODE = "chunkinfo";
var KEYNODE = "key";
var TOPICNODE = "topic";
var SCREENSNODE = "screens";
var GLOSSARYNODE = "glossary";
var ENTRYNODE = "entry";
var REMOTENODE = "remote";
var BREADCRUMBSNODE = "breadcrumbs";
var ITEMNODE = "item";
var CSHINFONODE = "csh-info";
var WINDOWLISTNODE = "windowlist";
var WINDOWNODE = "window";
var XCOORD = "x";
var YCOORD = "y";
var WIDTH = "width";
var HEIGHT = "height";
var OPTIONS = "options";

var REF = "ref";
var MASTERPROJECT = "MasterProject";
var CHILDID = "mergedchildid";

var LOADCOMPLETE = "rhloadcomplete";
var COOKIESPAGE = "access_cookies.htm";
var COOKIESPAGEID = "rhcookiereadwrite";

var DEFAULTURL = "defaulturl";
var MINWIDTH = "minwidth";
var MAXWIDTH = "maxwidth";
var MINHEIGHT = "minheight";
var MAXHEIGHT = "maxheight";
var BROWSERAGENT = "browseragent";
var FOLDER = "folder";
var VALUE = "value";
var DEFAULT = "default";
var RHMSWINNAME = "RHMSWINNAME";
var MAPNUM = "mapnum";
var MAPID = "mapid";
var TOPICURL = "topicurl";
var BCID = "#bc-";
var TOCCHILDIDPREFIX = "@";
var IDXLOADINGDIVID = "rhidxloadingdivid";
var GLOLOADINGDIVID = "rhgloloadingdivid";

var CSHMODE = "1";
var NONCSHMODE = "0";
var SHOWINCSHMODE = "CSH";
var SHOWINNONCSHMODE = "NONCSH";

var TRUESTR = "1";
var FALSESTR = "0";

var RHMAPID = "rhmapid";
var RHMAPNO = "rhmapno";
var RHWINDOW = "rhwnd";
var RHCSHMODE = "rhcsh";
var RHNEWWINDOW = "rhnewwnd";
var RHANDSEARCH = "rhandsearch";
var RHSEARCHSTR = "rhsearch";
var RHSYNSTR	= "rhsyns";
var RHHIGHLIGHT = "rhhl";
var RHSEARCHCOUNT = "rhsearchcount";
var RHHIGHLIGHTTEXTCOLOR = "rhhltxtcol";
var RHHIGHLIGHTBGCOLOR = "rhhlbgcol";

var DEFAULTXTHIGHLIGHTCOLOR = "#000000";
var DEFAULTBGHIGHLIGHTCOLOR = "#b2b4bf";

var WINLOCATION=0x01;		/*need location bar?*/
var WINMENUBAR=0x02;		/*need menubar?*/		
var WINRESIZABLE=0x04;	/*resizable window?*/
var WINTOOLBAR=0x08;		/*need toolbar?*/
var WINSTATUS=0x10;		/*need statusbar?*/
var WINSCROLLBARS=0x20;	/*need scrollbars?*/

var NAME = "name";
var URL = "url";
var SRC = "src";
var CHILDNAME = "childname";


var TREEITEMCLASS = "treeitem";
var LISTITEMCLASS = "listitem";
var UNCLICKABLECLASS = "unclickable";
var UNSELECTABLECLASS = "unselectable";
var NOLINKANCHORCLASS = "nolink";
var HLISTCLASS = "hlist";
var HANDCURSORCLASS = "handcursor";


var DATASRC = "data-src";
var DATAURL = "data-url";
var DATAPATH = "data-path";
var DATAROOTPATH = "data-rootpath";
var DATAITEMTYPE = "data-type";
var DATAHOVERIMGSRC = "data-hoverimgsrc";
var DATAIMGSRC = "data-imgsrc";
var DATASELIMGSRC = "data-selimgsrc";
var DATATERM = "data-term";
var DATAPANEID = "data-contentid";
var DATATABBUTTONID = "data-tabid";
var DATACLASS = "data-class";
var DATACLASSHOVER = "data-classhover";
var DATACLASSSEL = "data-classsel";
var DATASHOWIN = "data-showin";
var DATAPH = "data-placeholder";

var BOOKDELIM = ".";
var PAGEDELIM = "_";
var TABBUTTONID = "rhtabbuttonid";
var TOCID = "rhtocid";

var ITEMTYPEBOOKCLOSED	= 0;
var ITEMTYPEBOOKOPEN	= 1;
var ITEMTYPEPAGE		= 2;
var ITEMTYPEURL			= 3;
var ITEMTYPELOADING		= 4;
var ITEMTYPEBOOKCHILDS	= 5;
var ITEMTYPEICON		= 6;
var ITEMTYPEKW			= 7;
var ITEMTYPELINK		= 8;
var ITEMTYPESUBKW		= 9;
var ITEMTYPECATEGORY	= 10;
var ITEMTYPETERM		= 11;
var ITEMTYPEDEF			= 12;
var ITEMTYPEMAPNO		= 13;
var ITEMTYPEMAPID		= 14;

var SCR_NONE		= -1;
var SCR_INIT		= 0;
var SCR_CHILD_TOC	= 1;
var SCR_PARENT_TOC	= 2;
var SCR_CHILD_IDX	= 3;
var SCR_PARENT_IDX	= 4;
var SCR_CHILD_GLO	= 5;
var SCR_PARENT_GLO	= 6;
var SCR_PARENT_BC	= 7;
var SCR_PARENT_FTS	= 8;
var SCR_CHILD_FTS	= 9;
var SCR_PARENT_TOCSYNC = 10;
var SCR_CHILD_CSH	= 11;

var SEARCHPAGEWIDTHRATIO = 70;
var ECS_NONE			= 0;
var ECS_FTSREADY		= 1;
var ECS_SEARCHING		= 2;
var ECS_FATALERROR		= 3;
var ECS_CANCELED		= 4;
var ECS_SEARCHFAILED	= 5;
var ECS_FOUND			= 6;

var READ_REQ			= 1;
var SAVE_REQ			= 2;

var JS_TAGTOKEN			= 1;
var JS_TEXTTOKEN		= 3;

var LINK_NAME_MACRO = "{%LINK_NAME%}";
var LINK_URL_MACRO = "{%LINK_URL%}";
var SEARCH_SUMMARY_MACRO = "{%SEARCH_SUMMARY%}";
var SEARCH_URL_MACRO = "{%SEARCH_URL%}";
var ICON_MACRO = "{%ICON%}";