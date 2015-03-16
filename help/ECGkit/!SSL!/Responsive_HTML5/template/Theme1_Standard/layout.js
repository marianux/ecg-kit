/*************
	
	Custom setting names. Used to hook into RH settings mechanism.

************/
var lastVisitedTopic = "cust_last_visited_topic";
var sidebarCSHRunCheck = "cust_sidebar_csh_run_check";
/*************
	
	Global variables

************/
var showSidebarInCSHMode = false;
var desktopSidebar = null;
/*************
	
	Special actions for mobile devices.

************/
var iOS = ( navigator.userAgent.match(/(iPad|iPhone|iPod)/g) ? true : false );
if(iOS && gbBlockIOSScaling) {/* Only when scaling is disabled */
	window.onorientationchange = function()
	{
	   setTimeout("window.location.reload()", 500);
	}
}

var ua = navigator.userAgent;
if( ua.indexOf("Android") >= 0 )
{
	var Android = true;
	/* Refresh on orientation change to help problems with recalculating sidebar position */
	window.onorientationchange = function()
	{
	   setTimeout("window.location.reload()", 500);
	}
  
  /* Android 2.3 set scaling to avoid issues. */
  var androidversion = parseFloat(ua.slice(ua.indexOf("Android")+8)); 
  if (androidversion < 2.3)
  {
	var tmpMetaTagsList = document.getElementsByTagName('meta');
	var contentString = "user-scalable=no, initial-scale=1.0, maximum-scale=1.0, minimum-scale=1.0";
	for (i=0; i<tmpMetaTagsList.length; i++)
	{
		if (tmpMetaTagsList[i].name == "viewport")
			tmpMetaTagsList[i].content = contentString;
	}
	if(i == tmpMetaTagsList.length)
	{
		var metaTag = document.createElement('meta');
		metaTag.setAttribute("name", "viewport");
		metaTag.setAttribute("content", contentString);
		var headTags = document.getElementsByTagName('head');
		headTags[0].appendChild(metaTag);
	}
  }
}
/*************
	
	Initialize the layout

************/
  
$(function(){
	/* More global vars */
	$body = $("body");
	$sidebar = $("div.sidebar");
	$content = $("div.content");
	$contentholder = $("div.contentholder", $content);
	$sidebarButton = $("div.closebutton");
	
	/*Set localisation titles*/
	setLocalisation();

	/* Modify layout back to top link. RH changes the link on generation, set it correct now. */
	$("a", "div.content_top_link").attr("href", "#content_top_placeholder");
	
	/* Modify layout back to top link. RH changes the link on generation, set it correct now. */
	$topPlaceHolder = $("a", "div.content_top_link");
	$(window).load(function(){ contentTopToggle(); });
	$(window).resize(function(){ contentTopToggle(); });
	
	/* Use RH's event scheduling service. This service makes custom script work better with RH widgets. Adding scripts on DOM ready may interfere with widget loading. */
	addRhLoadCompleteEvent(initializeLayout);
	
	/* Set Flash objects to not overflow other elements */
	$("object").each(function(){
		$(this).attr("wmode", "transparent")
			   .append('<param name="wmode" value="transparent">')
			   .append('<embed wmode="transparent" height="'+$(this).attr('height')+'" width="'+$(this).attr('width')+'">');
	});
	
});
function initializeLayout() {/* This function makes changes to the layout, including widgets. Called through RH event handler. */
	
	/* Add hover event to default RH widgets (RH doesn't support it out of the box)*/
	hoverwidget($("div.home"));
	hoverwidget($("div.brsBack"));
	hoverwidget($("div.print"));
	hoverwidget($("div.brsNext"));
	
	readSetting(lastVisitedTopic, changeTopicLink, false, false);/* Check for the latest visited topic. */
	
	/* Make the sidebar resizable */
	desktopSidebar = new resizableSidebar($sidebar, "div.content", $body);
	desktopSidebar.setToggleElement($sidebarButton);
	desktopSidebar.toggleElement.click(function(){ desktopSidebar.toggle() });
	readSetting(desktopSidebar.lastSideBarWidthSetting, setLastSidebarWidth, false, false);/* Load sidebar width setting */
	readSetting(sidebarCSHRunCheck, setSidebarCSHRunCheck, false, false);/* Begin the CSH check for the sidebar. Only show/hide the sidebar in CSH mode the first time the help is opened. */
	$(window).resize(function() {
		if(desktopSidebar.sidebar.is(":visible")) {
			waitForFinalEvent(function(){
				desktopSidebar.reloadPane();
			}, 50, 'resize_content');
		} else {
			$sidebar.attr("style", "");
			$content.attr("style", "");
		}
	});
	
	setCheckBoxes();/* Change checkboxes into On/Off button */
	
	/* Set tab group hover */
	$("ul.wTabGroup li").hover(
							function(){tabButtonEvent($("ul.wTabGroup"))},
							function(){tabButtonEvent($("ul.wTabGroup"))})
						.click(function(){tabButtonEvent($("ul.wTabGroup"))});
						
	setTimeout(function(){ tabButtonEvent($("ul.wTabGroup"));/* Set images in tabs on sidebar */ }, 10);/* Set timeout because IE doesn't show the correct image directly. */
	
	var advancedSearch = new elementToggle($("div.cSearchSettingsToggle a"), $("div.wSearchResultSettings"));
	advancedSearch.setLinkState();
	advancedSearch.link.click(function(){ advancedSearch.toggle(); });
}
/*************
	
	Layout functions

************/

function contentTopToggle() {
	if($contentholder.prop('scrollHeight') > $contentholder.height() ) {
		$topPlaceHolder.show();
	} else {
		$topPlaceHolder.hide();
	}
}
var elementToggle = function(link, target) {
	this.link = link;
	this.target = target;
	this.linkVisibleClass = "opened";
	this.linkHiddenClass = "collapsed";
	this.state = function() { return (this.target.is(":visible")) ? true : false; };
	this.setLinkState = function (){
		if(this.state()) {
			this.link.removeClass(this.linkVisibleClass);
			this.link.addClass(this.linkHiddenClass);
		} else {
			this.link.addClass(this.linkVisibleClass);
			this.link.removeClass(this.linkHiddenClass);
		}
	}
	this.toggle = function() {
		if(this.state()) {
			this.target.hide();
			this.setLinkState();
		} else {
			this.target.show();
			this.setLinkState();
		}
	}
}
function hoverwidget(element) {
	 element.hover(
				function(){hoverWidgetButton($(this), 'hover')},
				function(){hoverWidgetButton($(this), 'unhover')})
			.click(function(){hoverWidgetButton($(this), 'hover')});
}
function hoverWidgetButton(element, action) {
	$hover = $("img.hover", element);
	$normal = $("img.normal", element);
	switch(action) {
		case 'hover':
				$hover.show();
				$normal.hide();
			break;
		case 'unhover':
				$normal.show();
				$hover.hide();
			break;s
	}
}
/*************
	
	Localisation

************/
function setLocalisation() {

	var localisationSetter = function(elementArray, string) {
		this.elements = elementArray;
		this.localisation = string;
		this.setLocalisation = function() {
			for(i=0;i<this.elements.length;i++) {
				var toSet = this.elements[i];
				toSet.attr("title", this.localisation).attr("alt", this.localisation);
				toSet.children().attr("title", this.localisation);
			}
		}
	}
	
	/* Set TOC title */
	var tocLocal = new localisationSetter(new Array($("div.bar_toc"), $("img", "#tocTabButton")), $("#localisation_toc").text());
	tocLocal.setLocalisation();
	var tocLocal = new localisationSetter(new Array($("div.bar_toc_active"), $("img", "#tocTabButton")), $("#localisation_toc").text());
	tocLocal.setLocalisation();
	/* Set IDX title */
	var idxLocal = new localisationSetter(new Array($("div.bar_index"), $("img", "#idxTabButton")), $("#localisation_idx").text());
	idxLocal.setLocalisation();
	var idxLocal = new localisationSetter(new Array($("div.bar_index_active"), $("img", "#idxTabButton")), $("#localisation_idx").text());
	idxLocal.setLocalisation();
	/* Set GLO title */
	var gloLocal = new localisationSetter(new Array($("div.bar_glossary"), $("img", "#gloTabButton")), $("#localisation_glo").text());
	gloLocal.setLocalisation();
	var gloLocal = new localisationSetter(new Array($("div.bar_glossary_active"), $("img", "#gloTabButton")), $("#localisation_glo").text());
	gloLocal.setLocalisation();
	/* Set FTS title */
	var ftsLocal = new localisationSetter(new Array($("div.bar_search")), $("#localisation_fts").text());
	ftsLocal.setLocalisation();
	var ftsLocal = new localisationSetter(new Array($("div.bar_search_active")), $("#localisation_fts").text());
	ftsLocal.setLocalisation();
	/* Set BACK title */
	var topicLocal = new localisationSetter(new Array($("div.bar_topic")), $("#localisation_back").text());
	topicLocal.setLocalisation();
}
/*************
	
	Sidebar

************/
var resizableSidebar = function(sidebar, reverseWidthElementName, baselineElement) {
	this.sidebar = sidebar;
	this.sidebarWidth = function() { return this.sidebar.width(); };
	this.setSidebarWidth = function(width) { this.sidebar.width(width); };
	this.reverse = reverseWidthElementName;
	this.reverseElement = $(reverseWidthElementName);
	this.reverseWidth = function() { return this.reverseElement.width(); };
	this.setReverseElementWidth = function(width) { this.reverseElement.width(width); };
	this.baseline = baselineElement;
	this.baselineWidth = function () { return this.baseline.width(); };
	this.maximumWidth = null
	this.minimumWidth = null;
	this.setMaximumWidth = function(maximum) { this.maximumWidth = maximum};
	this.resetMaximumWidth = function() { this.maximumWidth = this.baselineWidth()*0.75;};
	this.setMinimumWidth = function(minimum) { this.minimumWidth = minimum};
	this.resetMinimumWidth = function() { 
		if(this.toggleElement !== null) {
			this.minimumWidth = this.toggleElementWidth();
		} else {
			this.minimumWidth = this.baselineWidth()*0.02;
		}
	};
	this.handles = "e";
	
	this.toggleElement = null;
	this.toggleElementWidth = function () { return this.toggleElement.width(); };
	this.setToggleElement = function(element) { this.toggleElement = element; }
	this.toggleElementVisibleStyle = "buttonOpened";
	this.toggleElementHiddenStyle = "buttonClosed";
	this.setToggleOpeningState = function() { this.toggleElement.addClass(this.toggleElementHiddenStyle); this.toggleElement.removeClass(this.toggleElementVisibleStyle); };
	this.setToggleClosingState = function() { this.toggleElement.addClass(this.toggleElementVisibleStyle); this.toggleElement.removeClass(this.toggleElementHiddenStyle); };
	this.setToggleState = function() {
		if(this.sidebarWidth() > this.toggleElementWidth()) {
			this.setToggleClosingState();
		} else {
			this.setToggleOpeningState();
		}
	}
	this.toggle = function() {
		if(this.toggleElement === null) {
			return false;
		} else {
			if(this.sidebarWidth() > this.toggleElementWidth()) {
				this.setSidebarWidth(this.toggleElementWidth()-1);
			} else {
				this.setSidebarWidth(this.baselineWidth()*0.22);
			}
			this.setToggleState();
			this.reloadPane();
		}
	}
	
	this.lastSideBarWidthSetting = "cust_last_sidebar_width"+this.sidebar.reverseWidthElementName;
	this.saveSidebarWidth = function () { saveSetting(this.lastSideBarWidthSetting, this.sidebarWidth(), false); };
	
	this.reloadPane = function() {/* Function to reset pane sizes after page change or sidebar toggle */
		
		if(!this.sidebar.is(":visible")) {
			return false;/* Hidden, so do nothing. */
		}
		
		/* Set maximum width */
		this.resetMaximumWidth();
		this.resetMinimumWidth();
		this.sidebar.resizable({maxWidth: this.maximumWidth, minWidth: this.minimumWidth});
		
		if(this.baselineWidth() != (this.sidebarWidth() + this.reverseWidth())) {
			if(this.sidebarWidth() > this.baselineWidth()) {
				this.sidebar.width(this.maximumWidth);
			}
			
			this.reverseElement.width(this.baselineWidth() - this.sidebarWidth());
		}
		this.saveSidebarWidth();
		this.setToggleState();		
	}
	
	this.resetMaximumWidth();
	this.resetMinimumWidth();
	this.sidebar.resizable({alsoResizeReverseWidth: this.reverse, handles: this.handles, maxWidth: this.maximumWidth});
	
}
function setLastSidebarWidth(width, arg1, arg2) {
	if(width != "") {
		desktopSidebar.setSidebarWidth(width);
		desktopSidebar.reloadPane();
	}
}
function setSidebarCSHRunCheck(alreadyrunning, a, b) {
	return false;
	if(alreadyrunning == "") {
		saveSetting(sidebarCSHRunCheck, 'already_running', false);
		readSetting(RHCSHMODE, sidebarCSHCall);
	}
}
function sidebarCSHCall(cshmode, a, b) {/* Auto hide sidebar in CSH call */
	if(cshmode == CSHMODE) {
		if(showSidebarInCSHMode === false) {
			desktopSidebar.toggle();
		}
	}
}
/*************
	
	Tab buttons

************/
function tabButtonEvent(tabs) {
	$("li", tabs).each(function(){
		var img = $("img.normal", $(this));

		if($(this).attr("class") == "wTab") {
			hoverWidgetButton($(this), 'unhover');
		} else {
			hoverWidgetButton($(this), 'hover');
		}
	});
}
/*************
	
	Saving the last visited topic

************/
function setTopic() {
	var url = document.location.toString(),
		hashString = rh._.extractHashString(url);
	if(!hashString || -1 === hashString.search('ux=search')) {
		saveSetting(lastVisitedTopic, url, false);/* Save the URL of the topic. */
		changeTopicLink(url);
	}
}
function changeTopicLink(url, arg1, arg2) {/* Amend return link. Change default topic to last visited topic. */
	if(url != "") {
		$("a.back_to_last_topic").attr("href", url);
	}
}
/*************
	
	Turn checkboxes into sliders

************/
function setCheckBoxes() {
	/* 
		Credit for On/Off switch logic: http://blog.dbrgn.ch/2011/11/21/jquery-checkbox-on-off-toggleswitch/
		
		Method is shared under Creative Commons license Attribution-ShareAlike 3.0 Unported - http://creativecommons.org/licenses/by-sa/3.0/
		
	*/
	$('.switch').each(function() {
		var checkbox = $(this).children('input[type=checkbox]');
		var toggle = $(this).children('label.switch-toggle');
		if (checkbox.length) {
			
			checkbox.addClass('hidden');
			toggle.removeClass('hidden');
			
			if (checkbox[0].checked) {
				toggle.addClass('on');
				toggle.removeClass('off');
				toggle.text(toggle.attr('data-on'));
			} else { 
				toggle.addClass('off');
				toggle.removeClass('on');
				toggle.text(toggle.attr('data-off'));
			};
		}
	});
	$('.switch-toggle').click(function(){
		
		var checkbox = $(this).siblings('input[type=checkbox]')[0];
		var toggle = $(this); 
		
		if (checkbox.checked) {
			toggle.addClass('off');
			toggle.removeClass('on');
			toggle.text(toggle.attr('data-off'));
		} else {
			toggle.addClass('on');
			toggle.removeClass('off');
			toggle.text(toggle.attr('data-on'));
		}
	});
}
/*************
	
	Helper functions

************/
var waitForFinalEvent = (function () {
  var timers = {};
  return function (callback, ms, uniqueId) {
    if (!uniqueId) {
      //Ignore
    }
    if (timers[uniqueId]) {
      clearTimeout (timers[uniqueId]);
    }
    timers[uniqueId] = setTimeout(callback, ms);
  };
})();
/*************
	
	jQueryUI plugins

************/
$.ui.plugin.add("resizable", "alsoResizeReverseWidth", {

	start: function () {
		var that = $(this).data("ui-resizable"),
			o = that.options,
			_store = function (exp) {
				$(exp).each(function() {
					var el = $(this);
					el.data("ui-resizable-alsoresize-reverse", {
						width: parseInt(el.width(), 10), height: parseInt(el.height(), 10),
						left: parseInt(el.css("left"), 10), top: parseInt(el.css("top"), 10)
					});
				});
			};

		if (typeof(o.alsoResizeReverseWidth) === "object" && !o.alsoResizeReverseWidth.parentNode) {
			if (o.alsoResizeReverseWidth.length) { o.alsoResizeReverseWidth = o.alsoResizeReverseWidth[0]; _store(o.alsoResizeReverseWidth); }
			else { $.each(o.alsoResizeReverseWidth, function (exp) { _store(exp); }); }
		}else{
			_store(o.alsoResizeReverseWidth);
		}
	},

	resize: function (event, ui) {
		
		var that = $(this).data("ui-resizable"),
			o = that.options,
			os = that.originalSize,
			op = that.originalPosition,
			delta = {
				height: (that.size.height - os.height) || 0, width: (that.size.width - os.width) || 0,
				top: (that.position.top - op.top) || 0, left: (that.position.left - op.left) || 0
			},

			_alsoResizeReverseWidth = function (exp, c) {
				$(exp).each(function() {
					var el = $(this), start = $(this).data("ui-resizable-alsoresize-reverse"), style = {},
						/* Only resize the width of the element. */
						//css = c && c.length ? c : el.parents(ui.originalElement[0]).length ? ["width", "height"] : ["width", "height", "top", "left"];
						css = c && c.length ? c : el.parents(ui.originalElement[0]).length ? ["width"] : ["width", "left"];

					$.each(css, function (i, prop) {
						var sum = (start[prop]||0) - (delta[prop]||0);
						if (sum && sum >= 0) {
							style[prop] = sum || null;
						}
					});

					el.css(style);
				});
			};

		if (typeof(o.alsoResizeReverseWidth) === "object" && !o.alsoResizeReverseWidth.nodeType) {
			$.each(o.alsoResizeReverseWidth, function (exp, c) { _alsoResizeReverseWidth(exp, c); });
		}else{
			_alsoResizeReverseWidth(o.alsoResizeReverseWidth);
		}
	},

	stop: function () {
		$(this).removeData("resizable-alsoresize-reverse");
	}
});
