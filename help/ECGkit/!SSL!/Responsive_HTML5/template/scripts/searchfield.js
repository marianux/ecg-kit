var gSearchPageFilePath = "";
var gbGenerateForSP = 0;

gRootRelPath = ".";

addRhLoadCompleteEvent(initSearchFieldSubmit);

rh.model.subscribe(rh.consts('KEY_SEARCH_TERM'), function (searchTerm) {
	searchTerm = searchTerm || '';
	var searchBoxes = document.querySelectorAll(".wSearchField");
    for (i = searchBoxes.length-1; i > -1; i--) { 
		searchBoxes[i].value = searchTerm; 
	}
});
	
function searchHelp(e, searchBoxId, cshmode)
{
	if(e == null || e.keyCode == 13 || e.type == 'submit')
	{
		if(e != null) 
		{
			if(gbGenerateForSP || (e.type == 'submit' && cshmode == CSHMODE))
				preventEvent(e);
			e.target.blur();
		}
			
		var searchBox = document.getElementById(searchBoxId);
		var placeholderText = searchBox.getAttribute(DATAPH);
		if(searchBox == null || searchBox == 'undefined' || trimString(searchBox.value) == "" || (trimString(searchBox.value)==placeholderText && gbIE55 &&!gbIE10)) {
			return;
		}
		
		rh.model.publish(rh.consts('KEY_SEARCH_TERM'), searchBox.value);
		rh.model.publish(rh.consts('EVT_SEARCH_TERM'), true, {sync: true});
		return false;
	}
}
function initSearchFieldSubmit()
{
	if(gbIE5)
		readSetting(RHCSHMODE, callbackSearchFieldSubmit);
}
function callbackSearchFieldSubmit(cshmode)
{
	if(cshmode == CSHMODE && !gbPreviewMode)
	{
		var inputs = document.getElementsByTagName('input');
		for(var i=0; i<inputs.length; i++)
		{
			var searchAttr = inputs[i].getAttribute('data-search');
			if(searchAttr != null && searchAttr != 'undefined' && searchAttr == 'true')
			{
				var input = inputs[i];
				var id = input.getAttribute('id');
				patchInputForSubmit(input, function(){searchHelp(event, id, cshmode );});
			}
		}
	}
}