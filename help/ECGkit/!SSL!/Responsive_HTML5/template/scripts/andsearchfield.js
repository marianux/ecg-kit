var gANDSearchEnabled = 0;
var gANDSearchControlId = "andsearch";

initANDSearchControl();
function initANDSearchControl()
{
	addRhLoadCompleteEvent(loadANDSearchSetting);
}
function loadANDSearchSetting()
{
	readSetting(RHANDSEARCH, callbackAndSearchCookieRead);
}
function callbackAndSearchCookieRead(andSearchFlag)
{
	if(andSearchFlag == TRUESTR)
		gANDSearchEnabled = 1;
	else if(andSearchFlag == FALSESTR)
		gANDSearchEnabled = 0;
	var andSearchElem = document.getElementById(gANDSearchControlId);
	if(andSearchElem)
	{
		if(gANDSearchEnabled)
			andSearchElem.checked = true;
		else
			andSearchElem.checked = false;
	}
}
function onToggleANDSearch()
{
	if(gANDSearchEnabled)
	{
		gANDSearchEnabled = 0;
		saveSetting(RHANDSEARCH, FALSESTR, true);
	}
	else
	{
		gANDSearchEnabled = 1;
		saveSetting(RHANDSEARCH, TRUESTR, true);
	}
}