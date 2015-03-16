var gbCookieSupported=null;
var gbCookieSupportedWithoutPath=null;
var gbLocalStorageSupported=null;

var gHost = null;
var gHostPath = "/";
var gIFrameElem = null;

function saveSetting(name, value, bPersistent)
{
	if(isCookieFullySupported() || gHost == null)
		setCookie(name, value, bPersistent);
	else if(isLocalDBSupported())
		setInLocalDB(name, value, bPersistent);
	else if(isCookieSupportedWithoutPath())
		setThroughIFrame(name, value, bPersistent);
}
function readSetting(name, oCallbackFunc, arg1, arg2)
{
	var val="";
	if(isCookieFullySupported() || gHost == null)
	{
		val = getCookie(name);
		if(oCallbackFunc)
			oCallbackFunc(val, arg1, arg2);
	}
	else if(isLocalDBSupported())
	{
		val = getFromLocalDB(name);
		if(oCallbackFunc)
			oCallbackFunc(val, arg1, arg2);
	}
	else if(isCookieSupportedWithoutPath())
		getThroughIFrame(name, oCallbackFunc, arg1, arg2);
	else
		return false;
		
	return true;
}
function isCookieFullySupported()
{
	if(gbCookieSupported)
		return gbCookieSupported;
	document.cookie = "testcookieone;domain=" + gHost + ";path="+gHostPath;
	if(document.cookie.indexOf("testcookieone")!=-1)
		gbCookieSupported = true;
	else
		gbCookieSupported = false;
	return gbCookieSupported;
}
function isCookieSupportedWithoutPath()
{
	if(gbCookieSupportedWithoutPath)
		return gbCookieSupportedWithoutPath;
	document.cookie = "testcookietwo";
	if(document.cookie.indexOf("testcookietwo")!=-1)
		gbCookieSupportedWithoutPath = true;
	else
		gbCookieSupportedWithoutPath = false;
	return gbCookieSupportedWithoutPath;
}
function isLocalDBSupported()
{
	if(gbLocalStorageSupported === null) {
		if(localStorage != undefined) {
			try {
				localStorage.setItem('dummyTestForLocalStorage', 1);
				localStorage.removeItem('dummyTestForLocalStorage');
			} catch(e) {
				gbLocalStorageSupported = false;
			}

			if(gbLocalStorageSupported === null) {
				gbLocalStorageSupported = true;
			}
		} else {
			gbLocalStorageSupported = false;
		}		
	}
	return gbLocalStorageSupported;
}
function initSettings(commonRootRelPath)
{
	if(commonRootRelPath == null || commonRootRelPath == "")
		return;
	var folderAbsPath = _getFullPath(_getPath(document.location.href), commonRootRelPath + "/");
	if(_isHTTPUrl(folderAbsPath))
	{
		gHost = _getHostNameFromURL(folderAbsPath);
		gHostPath = _getPathFromURL(folderAbsPath);
	}
	else
	{
		gHost = "";
		gHostPath = folderAbsPath;
	}
}
function setCookie(name, value, bPersistent)
{
	var expires = ";";
	if(bPersistent)
		expires = ";expires=Sunday, March 15, 2037 9:30:01 AM;";
	if(gHost!=null && gHostPath!=null)
		document.cookie = name + "=" + value + ";domain=" + gHost + ";path=" + gHostPath + expires;
	else
		document.cookie = name + "=" + value;
}
function getCookie(name)
{

	var val = "";
	var delimiter = ";";
	if (document.cookie == '')
		return '';
	else
	{
		var namePos, valueStartPos, valueEndPos;
		var rawCookie = document.cookie;
		while(rawCookie.length > 0) {
			namePos = rawCookie.indexOf(name);
			if(namePos == -1) {
				break;
			}
			
			var equalPos = namePos + name.length;
			if(rawCookie.charAt(equalPos) == '=')
			{
				valueStartPos = equalPos + 1;
				valueEndPos = rawCookie.indexOf(delimiter, valueStartPos);
				if (valueEndPos == -1) valueEndPos = rawCookie.length;
				val = rawCookie.substring(valueStartPos, valueEndPos);
				break;
			}
			rawCookie = rawCookie.substring(equalPos);
		}
	}
	
	if(val == null || val == 'undefined')
		val = '';
	
	return val;
}
function setInLocalDB(name, value, bPersistent)
{
	try {
		if(bPersistent)
			localStorage.setItem(name + gHost + gHostPath, value);
		else
			sessionStorage.setItem(name + gHost + gHostPath, value);
	} catch(e) {}
}
function getFromLocalDB(name)
{
	var val = '';
	try {
		val = sessionStorage.getItem(name + gHost + gHostPath);
		if(val == null)
			val = localStorage.getItem(name + gHost + gHostPath);
	} catch (e) {}
			
	if(val == null || val == undefined)
		val = '';
	return val;
}

var cookieRequestQ = new MhQueue();
var gbIFrameLoaded = false;
var gbIFrameLoading = false;
function cookieSaveRequesObj(reqType, name, value, bPersistent)
{
	this.reqType = reqType;
	this.name = name;
	this.value = value;
	this.bPersistent = bPersistent;
}
function cookieReadRequesObj(reqType, name, oCallbackFunc, arg1, arg2)
{
	this.reqType = reqType;
	this.name = name;
	this.oCallbackFunc = oCallbackFunc;
	this.arg1 = arg1;
	this.arg2 = arg2;
}
function insertIFrame()
{
	gbIFrameLoading = true;
	gIFrameElem = document.createElement('iframe');
	gIFrameElem.setAttribute("src", gHostPath+COOKIESPAGE);
	gIFrameElem.setAttribute("name", COOKIESPAGEID);
	gIFrameElem.setAttribute("height", "0px");
	gIFrameElem.setAttribute("width", "0px");
	gIFrameElem.style.display = "none";
	var bodyElem = document.getElementsByTagName('body')[0];
	bodyElem.appendChild(gIFrameElem);
	if( gIFrameElem.addEventListener )
		gIFrameElem.addEventListener('load', performRequest, false);
	else if( gIFrameElem.attachEvent )
		gIFrameElem.attachEvent('onload', performRequest, false);	
}
function setThroughIFrame(name, value, bPersistent)
{
	var objSave = new cookieSaveRequesObj(SAVE_REQ, name, value, bPersistent);
	cookieRequestQ.enqueue(objSave);
	
	if(!gbIFrameLoaded)
	{
		if(gbIFrameLoading)
			return;
		else
			insertIFrame();
	}
	else
		performRequest();
}

function getThroughIFrame(name, oCallbackFunc, arg1, arg2)
{
	var objRead = new cookieReadRequesObj(READ_REQ, name, oCallbackFunc, arg1, arg2)
	cookieRequestQ.enqueue(objRead);
	
	if(!gbIFrameLoaded)
	{
		if(gbIFrameLoading)
			return;
		else
			insertIFrame();
	}
	else
		performRequest();
}
function performRequest()
{
	gbIFrameLoaded = true;
	gbIFrameLoading = false;
	if(cookieRequestQ.isEmpty())
		return;
	var obj = cookieRequestQ.dequeue();
	if(obj.reqType == SAVE_REQ)
		gIFrameElem.contentWindow.setIFrameCookie(obj.name, obj.value, obj.bPersistent, gHost, gHostPath);
	else if(obj.reqType == READ_REQ)
	{
		var val = gIFrameElem.contentWindow.getIFrameCookie(obj.name);
		if(obj.oCallbackFunc)
			obj.oCallbackFunc(val, obj.arg1, obj.arg2);
	}
	performRequest();
}