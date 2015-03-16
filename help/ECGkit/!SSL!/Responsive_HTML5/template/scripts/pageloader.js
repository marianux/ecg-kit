initPage();
function initPage()
{
	if( window.addEventListener )
	{
		window.addEventListener('load', onPageLoad, false);
	}
	else if( window.attachEvent )
	{
		window.attachEvent('onload', onPageLoad, false);
	}
	if( window.addEventListener )
	{
		window.addEventListener('unload', onPageUnload, false);
	}
	else if( window.attachEvent )
	{
		window.attachEvent('onunload', onPageUnload, false);
	}

}
function onPageUnload()
{
	//This function is deliberately set empty. 
	//Overriding onunload results in clearing of back-forward cache(also called bfcache)
}

function onPageLoad()
{
	if(gbBlockIOSScaling)
		blockIOSScaling();
	if(gbPreviewMode)
	{
		initSettings(gRootRelPath);
		fireRhLoadCompleteEvent();
	}
	else
		initAndLoadParentData(null, SCR_NONE);
}