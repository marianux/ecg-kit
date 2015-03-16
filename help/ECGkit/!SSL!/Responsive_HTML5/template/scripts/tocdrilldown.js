var gTopicId = "";
function DDList(rootRelPath, dataFolder, rootFile, commonRootRelPath)
{
	this.errorMsg = "";
	this.closedBookClass = "";
	this.openBookClass = "";
	this.pageClass = "";
	this.urlClass = "";
	this.closedBookClassHover = "";
	this.openBookClassHover = "";
	this.pageClassHover = "";
	this.urlClassHover = "";
	this.closedBookClassClick = "";
	this.openBookClassClick = "";
	this.pageClassClick = "";
	this.urlClassClick = "";
	this.closedBookClassSelected = "";
	this.openBookClassSelected = "";
	this.pageClassSelected = "";
	this.urlClassSelected = "";
	this.pageHtml = "";
	this.bookClosedHtml = "";
	this.bookOpenHtml = "";
	this.urlHtml = "";
	
	this.iconClass = new Object();
	this.iconClass[ITEMTYPEBOOKCLOSED] = "";
	this.iconClass[ITEMTYPEBOOKOPEN] = "";
	this.iconClass[ITEMTYPEPAGE] = "";
	this.iconClass[ITEMTYPEURL] = "";
	
	this.iconStyle = new Object();
	this.iconStyle[ITEMTYPEBOOKCLOSED] = "";
	this.iconStyle[ITEMTYPEBOOKOPEN] = "";
	this.iconStyle[ITEMTYPEPAGE] = "";
	this.iconStyle[ITEMTYPEURL] = "";
	
	this.iconSrc = new Object();
	this.iconSrc[ITEMTYPEBOOKCLOSED] = "";
	this.iconSrc[ITEMTYPEBOOKOPEN] = "";
	this.iconSrc[ITEMTYPEPAGE] = "";
	this.iconSrc[ITEMTYPEURL] = "";
	
	this.iconHoverSrc = new Object();
	this.iconHoverSrc[ITEMTYPEBOOKCLOSED] = "";
	this.iconHoverSrc[ITEMTYPEBOOKOPEN] = "";
	this.iconHoverSrc[ITEMTYPEPAGE] = "";
	this.iconHoverSrc[ITEMTYPEURL] = "";
	
	this.iconSelSrc = new Object();
	this.iconSelSrc[ITEMTYPEBOOKCLOSED] = "";
	this.iconSelSrc[ITEMTYPEBOOKOPEN] = "";
	this.iconSelSrc[ITEMTYPEPAGE] = "";
	this.iconSelSrc[ITEMTYPEURL] = "";	
	
	this.iconHtml = new Object();
	this.iconHtml[ITEMTYPEBOOKCLOSED] = "";
	this.iconHtml[ITEMTYPEBOOKOPEN] = "";
	this.iconHtml[ITEMTYPEPAGE] = "";
	this.iconHtml[ITEMTYPEURL] = "";		
	
	this.upButtonId = "";

	this.rootFile = rootFile;
	this.rootRelPath = rootRelPath;
	this.commonRootRelPath = commonRootRelPath;	
	this.dataFolder = dataFolder;
	this.rootHtmlNode = null;
	
	this.selectedListNode = null;
	this.hoveredListNode = null;
	this.syncList = true;
	
	this.urlId = ""
	this.idParts = null;
	this.idPartsNextIndex = 0;
	this.curParentId = "";
	this.nextNodeIdToBeSynced = null;
	this.doSyncNeeded = true;
	
	this.loadingIconClass = "loadingicon";
	this.loadingIconHtml = "";
	this.loadingTextClass = "loadingtext";
	this.loadingText = "";
	this.isSyncingRoot = true;
	this.childProjOrder = 0;
	
	this.isBookChildScreen = false;
	
	this.loadStack = new MhStack();
	this.parentBookSrcStack = new MhStack();
	
	DDList.prototype.init = function()
	{
		this.rootHtmlNode = document.getElementById(this.rootId);
		if(this.rootHtmlNode.attachEvent)
		{
			this.rootHtmlNode.attachEvent("onkeydown", function(){onKeyPress(event);});
			this.rootHtmlNode.attachEvent("onblur", function(){onListBlur();});
			this.rootHtmlNode.attachEvent("onfocus", function(){onListFocus();});
		}
		else
		{
			this.rootHtmlNode.setAttribute("onkeydown" , "onKeyPress(event)");	
			this.rootHtmlNode.setAttribute("onblur" , "onListBlur()");	
			this.rootHtmlNode.setAttribute("onfocus" , "onListFocus()");	
		}
	}
	DDList.prototype.load = function()
	{
		this.insertLoadingMsg();
		var parentId = this.rootHtmlNode.getAttribute('id');
		this.curParentId = parentId;
		var objContext = new listContext(0,0,0);
		var parentBkInfo = new parentBookInfo(this.rootRelPath, this.rootFile);
		this.parentBookSrcStack.push(parentBkInfo);
		this.updateUpButton(objContext);
		var srcPath = this.rootRelPath + "/" + this.dataFolder + "/" + this.rootFile;
		xmlJsReader.loadFile(srcPath, callbackCreateList, objContext);
	}
	DDList.prototype.createChildList = function(xmlDoc, objContext)
	{
		var dataXmlNode = xmlDoc.getElementsByTagName(DATANODE)[0];
		var parentUrl = "";	
		var parentName = dataXmlNode.getAttribute(NAME);
		if(parentName != null && parentName != 'undefined' && parentName != "")
		{
			var parentUrl = dataXmlNode.getAttribute(URL);
			if(parentUrl != "" && parentUrl != null && parentUrl != 'undefined')
				this.insertListItem(this.curParentId, parentName, parentUrl, objContext.bookCount, objContext.pageCount, objContext.childProjOrder, ITEMTYPEBOOKOPEN);
		}		
		this.insertListItems(dataXmlNode, 0, objContext);
	}
	DDList.prototype.insertListItems = function(dataXmlNode, index, objContext)
	{
		var bookCount=objContext.bookCount;
		var pageCount=objContext.pageCount;
		var parentId = this.curParentId;
		var origRootRelPath = objContext.rootRelPath;
		var origCommonRootRelPath = objContext.commonRootRelPath;
		var rootRelPath = origRootRelPath;
		var commonRootRelPath = origCommonRootRelPath;
		if(commonRootRelPath == null)
			commonRootRelPath = this.commonRootRelPath;
		if(rootRelPath == null)
			rootRelPath = this.rootRelPath;	
		var len = dataXmlNode.childNodes.length;

		for(var i = index; i < len; i++)
		{
			var xmlNode = dataXmlNode.childNodes[i];
			var name = xmlNode.getAttribute(NAME);
			var url = xmlNode.getAttribute(URL);
			if(url != null && url != 'undefined' && xmlNode.tagName != URLNODE && !_isRemoteUrl(url))
				url = rootRelPath + "/" + url;
			if(xmlNode.tagName == BOOKNODE)
			{
				bookCount++;
				pageCount = 0;
				objContext.bookCount = bookCount;
				objContext.pageCount = pageCount;
				var src = xmlNode.getAttribute(SRC);
				this.insertListItem(parentId, name, url, bookCount, pageCount, objContext.childProjOrder, ITEMTYPEBOOKCLOSED, src, origRootRelPath, origCommonRootRelPath);
			}
			else if(xmlNode.tagName == PAGENODE)
			{
				pageCount++;
				objContext.pageCount = pageCount;
				this.insertListItem(parentId, name, url, bookCount, pageCount, objContext.childProjOrder, ITEMTYPEPAGE);
			}
			else if(xmlNode.tagName == URLNODE)
			{
				pageCount++;
				objContext.pageCount = pageCount;
				this.insertListItem(parentId, name, url, bookCount, pageCount, objContext.childProjOrder, ITEMTYPEURL);
			}
			else if(xmlNode.tagName == PROJNODE)
			{
				++this.childProjOrder;
				objContext.childProjOrder = this.childProjOrder;
				var dInfo = new dataInfo(dataXmlNode, i+1, origRootRelPath, origCommonRootRelPath, bookCount, pageCount);
				this.loadStack.push(dInfo);
				objContext.bookCount = 0;
				objContext.pageCount = 0;
				var ref = xmlNode.getAttribute(REF);
				this.insertChildProjList(ref, objContext);
				break;
			}
		}
		if(len == i)
			this.loadFromStack(objContext);

	}
	DDList.prototype.loadFromStack = function(objContext)
	{
		objContext.childProjOrder = 0;
		if(this.loadStack.isEmpty() != true)
		{
			var dInfo =	this.loadStack.pop();
			objContext.rootRelPath = dInfo.rootRelPath;
			objContext.commonRootRelPath = dInfo.commonRootRelPath;
			objContext.bookCount = dInfo.bookCount;
			objContext.pageCount = dInfo.pageCount;
			this.insertListItems(dInfo.dataXmlNode, dInfo.index, objContext);
		}
		else
		{	
			this.childProjOrder = 0;		
			this.removeLoadingMsg();
			this.hoveredListNode = this.getFirstListNode();
			//this.focusHoveredNode();
			if(this.syncList == true)
				loadParentDataForSyncing(gCommonRootRelPath, SCR_PARENT_TOCSYNC);
		}
	}
	DDList.prototype.insertListItem = function(parentId, name, url, bookCount, pageCount, childProjOrder, itemType, src, rootRelPath, commonRootRelPath)
	{
		var listNode = document.createElement("div");
		listNode.setAttribute('class', LISTITEMCLASS);
		var classNormal = "";
		var classHover = "";
		var classClick = "";
		var html = "";
		var dot = "";
		var expandableNode = false;
		if(this.rootHtmlNode.getAttribute('id') != parentId)
			dot = BOOKDELIM;
		var strChildProjOrder = "";
		if(childProjOrder != 0)
			strChildProjOrder = "C" + childProjOrder;
		if(ITEMTYPEBOOKOPEN == itemType)
		{
			if(url != '' && url != null)
				listNode.setAttribute(DATAURL, url);
			listNode.setAttribute('id', parentId + strChildProjOrder);
			classNormal = this.openBookClass;
			if(url == null)
				classNormal = classNormal + " " + UNCLICKABLECLASS;
			classHover = this.openBookClassHover;
			classClick = this.openBookClassClick;
			html = this.bookOpenHtml.replace(LINK_NAME_MACRO, name);
			expandableNode = false;
		}
		else if(ITEMTYPEBOOKCLOSED == itemType)
		{
			if(src != null && src != '')
				listNode.setAttribute(DATASRC, src);
			if(url != '' && url != null)
				listNode.setAttribute(DATAURL, url);
			if(rootRelPath != null)
				listNode.setAttribute(DATAPATH, rootRelPath);
			if(commonRootRelPath != null)
				listNode.setAttribute(DATAROOTPATH, commonRootRelPath);				
			listNode.setAttribute('id', parentId + dot + bookCount + strChildProjOrder);
			classNormal = this.closedBookClass;
			if(url == null)
				classNormal = classNormal + " " + UNCLICKABLECLASS;
			classHover = this.closedBookClassHover;
			classClick = this.closedBookClassClick;
			html = this.bookClosedHtml;
			expandableNode = true;
		}
		else if(ITEMTYPEPAGE == itemType)
		{
			listNode.setAttribute('id', parentId + dot + bookCount + PAGEDELIM + pageCount+ strChildProjOrder);
			classNormal = this.pageClass;
			classHover = this.pageClassHover;
			classClick = this.pageClassClick;
			html = this.pageHtml;
			expandableNode = false;
		}
		else if(ITEMTYPEURL == itemType)
		{
			classNormal = this.urlClass;
			classHover = this.urlClassHover;
			classClick = this.urlClassClick;	
			html = this.urlHtml;
			expandableNode = false;
		}
		this.rootHtmlNode.appendChild(listNode);
		this.insertChildHtmlNode(listNode, name, itemType, html, classNormal, classHover, classClick, url, expandableNode);
	}
	DDList.prototype.insertChildHtmlNode = function(listNode, name, itemType, html, classNormal, classHover, classClick, url, expandableNode)
	{
		var bAddAnchor = false;
		if(url != null && url != "")
			bAddAnchor = true;
		html = html.replace(LINK_NAME_MACRO, name);
		var iconHtml = this.getIconHtml(itemType);
		html = html.replace(ICON_MACRO, iconHtml);
		var htmlNode = document.createElement("div");
		htmlNode.className = classNormal + " " + UNSELECTABLECLASS;
		htmlNode.setAttribute("title", name);
		this.setNodeItemType(htmlNode, itemType);
		htmlNode.innerHTML = html;
		if(bAddAnchor)
		{
			var anchorNode = document.createElement("a");
			anchorNode.className = NOLINKANCHORCLASS;
			anchorNode.appendChild(htmlNode);
			listNode.insertBefore(anchorNode, listNode.firstChild);
		}
		else
			listNode.insertBefore(htmlNode, listNode.firstChild);
		var urlWithId = this.addEventsToNode(htmlNode, classNormal, classHover, classClick, url, expandableNode);
		if(bAddAnchor)
		{
			if(urlWithId == "")
				url = "#";
			else if(getUrlBookmark(url) != "")
				url = urlWithId;
			else
				tocid = getUrlParameter(TOCID, urlWithId);
			anchorNode.setAttribute("href", url);	
		}
	}	
	DDList.prototype.insertChildProjList = function(ref, objContext)
	{
		var commonRootRelPath = objContext.commonRootRelPath;
		if(commonRootRelPath == null)
			commonRootRelPath = this.commonRootRelPath;
		var strChildProjPath = commonRootRelPath + "/" + ref;
		loadScreenData(strChildProjPath, SCR_CHILD_TOC, objContext);	
	}	
	DDList.prototype.addEventsToNode = function(htmlNode, classNormal, classHover, classClick, url, expandableNode)
	{
		if(htmlNode.attachEvent)
		{
			if(isTouchDevice())
			{
				htmlNode.attachEvent('ontouchstart', function(){onNodeHover(htmlNode, classHover);});
				htmlNode.attachEvent('ontouchend', function(){onNodeHoverOut(htmlNode, classHover);});
				htmlNode.attachEvent('ontouchmove', function(){onNodeHoverOut(htmlNode, classHover);});
			}
			else
			{
				htmlNode.attachEvent('onmouseout', function(){onNodeHoverOut(htmlNode, classNormal);});
				htmlNode.attachEvent('onmouseover', function(){onNodeHover(htmlNode, classHover);});
			}
		}
		else
		{
			if(isTouchDevice())
			{
				htmlNode.setAttribute("ontouchstart", "onNodeHover(this,'" + classHover + "')");
				htmlNode.setAttribute("ontouchend", "onNodeHoverOut(this,'" + classNormal + "')");
				htmlNode.setAttribute("ontouchmove", "onNodeHoverOut(this,'" + classNormal + "')");
			}
			else
			{
				htmlNode.setAttribute("onmouseout", "onNodeHoverOut(this,'" + classNormal + "')");
				htmlNode.setAttribute("onmouseover", "onNodeHover(this,'" + classHover + "')");
			}
		}
		
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		var id = listNode.getAttribute('id');
		var urlWithId = "";
		if(expandableNode == false)
		{
			if(url == null || url =='')
				return "";
			if(id != null && id != '')
			{
				var urlid = "";
				if(url.indexOf("?") != -1)
					urlid = '&' + TOCID + '=' + id;
				else
					urlid = '?' + TOCID + '=' + id;
				var bookmark = getUrlBookmark(url);
				if(bookmark != "")
					url = getUrlWithoutBookmark(url);
				urlWithId = url + urlid + bookmark;
			}
			else
				urlWithId = url;
			/*
			if(htmlNode.attachEvent)
				htmlNode.attachEvent('onclick', function(){onNodeClick(htmlNode, classClick, classNormal, urlWithId)});
			else
				htmlNode.setAttribute("onclick", "onNodeClick(this, '" + classClick + "', '" + classNormal + "', '" + urlWithId + "')");
			*/
		}
		else
		{
			if(htmlNode.attachEvent)
				htmlNode.attachEvent('onclick', function(){onBookClick(htmlNode)});
			else
				htmlNode.setAttribute("onclick", "onBookClick(this)");		
		
		}
		return urlWithId;
	}
	DDList.prototype.clickUp = function()
	{
		var upHtmlNode = this.getUpHtmlNode();
		if(upHtmlNode == null || upHtmlNode == 'undefined')
			return;
		var newParentId = this.getParentsParentId(this.curParentId);
		var url = this.prepareUrl(newParentId);
		document.location.href = url;
	}
	DDList.prototype.getParentsParentId  = function(parentId)
	{
		var newParentId = parentId;
		newParentId = newParentId.substring(0, newParentId.lastIndexOf("."));
		if(newParentId == "")
			newParentId = this.rootHtmlNode.getAttribute('id');
		return newParentId;
	}
	DDList.prototype.isRoot = function()
	{
		var id = this.curParentId;
		if(id == this.rootHtmlNode.getAttribute('id'))
			return true;
		else
			return false;
	}
	DDList.prototype.updateUpButton = function()
	{
		var upHtmlNode = this.getUpHtmlNode();
		if(upHtmlNode == null || upHtmlNode == 'undefined')
			return;
		if(this.isRoot())
		{
			this.isBookChildScreen = false;
			style = "none";
		}
		else
		{
			this.isBookChildScreen = true;
			style = "block";
		}
	
		upHtmlNode.style.display = style;
	}
	DDList.prototype.getUpHtmlNode = function()
	{
		var upButtonNode = document.getElementById(this.upButtonId)
		return upButtonNode;
	}
	DDList.prototype.sync = function(childPrefix, childOrder)
	{
		if(this.doSyncNeeded != true)
			return;
		if(this.isRoot())
		{
			this.urlId = getUrlParameter(TOCID, document.location.href);
			if(this.urlId == '')
			{
				if(gTopicId != null && gTopicId != 'undefined' && gTopicId != '')
				{
					var dataId = this.rootHtmlNode.getAttribute('id');
					var pos = gTopicId.indexOf(".");
					if(childOrder != "" )
						childOrder = "C" + childOrder;
					else
						childOrder = "";
					if(pos != -1)
					{
						var firstPart = gTopicId.substring(0, pos);
						var secondPart = gTopicId.substring(pos, gTopicId.length);
						this.urlId = dataId + childPrefix + firstPart + childOrder + secondPart;
					}
					else
						this.urlId = dataId + childPrefix + gTopicId + childOrder;
				}
				else
				{
					this.doSyncNeeded = false;
					return;
				}
			}
			this.idParts = this.urlId.split(BOOKDELIM);;
			this.nextNodeIdToBeSynced = this.idParts[this.idPartsNextIndex];
		}
		var listNode = this.getListNodeById(this.nextNodeIdToBeSynced);
		++this.idPartsNextIndex;
		if(this.idPartsNextIndex != this.idParts.length)
		{
			this.loadBook(listNode);
			this.nextNodeIdToBeSynced = this.nextNodeIdToBeSynced + BOOKDELIM + this.idParts[this.idPartsNextIndex];
		}
		else
		{
			this.doSyncNeeded = false;
			if(this.isBookListNode(listNode))
				this.loadBook(listNode);
			else
				this.setSelectedListNode(listNode);
		}
	}	
	DDList.prototype.pressKey = function(e)
	{
		if (e.preventDefault)
            e.preventDefault();
		var kCode = 0;
		if(e.keyCode)
			kCode = e.keyCode;
		else
			kCode = e.which;
		var listNode = null;
		var htmlNode = null;
		var event = "";
		if(kCode == 38)
		{
			listNode = this.getPreviousListItem(this.hoveredListNode);
			htmlNode = this.getHtmlNodeFromListNode(listNode);
			event = "mouseover";
		}
		else if(kCode == 40)
		{
			listNode = this.getNextListItem(this.hoveredListNode);
			htmlNode = this.getHtmlNodeFromListNode(listNode);
			event = "mouseover";                
		}
		else if(kCode == 13 || kCode == 32)
		{
			listNode = this.hoveredListNode;
			htmlNode = this.getHtmlNodeFromListNode(listNode);
			event = "click";
		}
		else if(kCode == 39)
		{
			if(this.isBookListNode(this.hoveredListNode))
			{
				listNode = this.hoveredListNode;
				htmlNode = this.getHtmlNodeFromListNode(listNode);
				event = "click";
			}
		}
		else if(kCode == 37)
		{
			if(this.isBookChildScreen)
			{
				htmlNode = this.getUpHtmlNode();
				event = "click";
			}
		}
		if(htmlNode != null)
			fireEvent(htmlNode, event);
	}	
	DDList.prototype.getListNodeById = function(id)
	{
		return document.getElementById(id);
	}	
	DDList.prototype.getHtmlNodeFromListNode = function(listNode)
	{
		if(listNode == null || listNode.nodeType != JS_TAGTOKEN)
			return null;
		var childHtmlNodes = null;
		var childs = listNode.childNodes;
		var len = childs.length;
		var i = 0;
		var anchorElem = null;
		for(i=0; i<len; i++)
		{
			anchorElem = childs[i];
			if(anchorElem.nodeType == JS_TAGTOKEN)
				break;	
		}
		if(anchorElem != null && anchorElem.nodeName == "A")
			childHtmlNodes = anchorElem.getElementsByTagName("div");
		else
			childHtmlNodes = listNode.getElementsByTagName("div");
		var htmlNode = null;
		htmlNode = childHtmlNodes[0];
		return htmlNode;		
	}
	DDList.prototype.getListNodeFromHtmlNode = function(htmlNode)
	{
		if(htmlNode != null && htmlNode != 'undefined')
		{
			var pNode = htmlNode.parentNode;
			if(pNode != null && pNode.nodeName == "A")
				pNode = pNode.parentNode;			
			if(pNode == this.rootHtmlNode)
				return null;
			else
				return pNode;
		}
		else
			return null;
	}	
	DDList.prototype.getIconHtmlNodeFromListNode = function(listNode)
	{
		if(listNode == null)
			return null;
		var htmlNode = this.getHtmlNodeFromListNode(listNode)
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		return iconHtmlNode;
	}
	DDList.prototype.getHtmlChildNode = function(node, tag, type)
	{
		if(tag == "" || tag == 'undefined')
			return null;
		var nodeChilds = node.getElementsByTagName(tag);
		var len = nodeChilds.length;
		var i=0;
		for(i=0; i<len; i++)
		{
			if(this.isNodeItemTypeThis(nodeChilds[i],type))
			 return nodeChilds[i];		
		}
		
	}
	DDList.prototype.getFirstListNode = function()
	{
		var listNodes = this.rootHtmlNode.getElementsByTagName("div");
		var i=0;
		var len = listNodes.length;
		var htmlNode = null;
		var itemType = -1;
		for(i=0; i<len; i++)
		{
			htmlNode = this.getHtmlNodeFromListNode(listNodes[i]);
			if(htmlNode != null)
			{
				itemType = this.getNodeItemType(htmlNode);
				if(itemType == ITEMTYPEBOOKOPEN ||
					itemType == ITEMTYPEBOOKCLOSED ||
					itemType == ITEMTYPEPAGE ||
					itemType == ITEMTYPEURL)
					return listNodes[i];
			}
		}
		return null;
	}	
	DDList.prototype.getNextListItem = function(listNode)
	{
		if(listNode == null || listNode == 'undefined')
			return null;
		var nextSibling = listNode.nextSibling;
		if(nextSibling == null || nextSibling.tagName == null || nextSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return nextSibling	
	}
	DDList.prototype.getPreviousListItem = function(listNode)
	{
		if(listNode == null || listNode == 'undefined')
			return null;
		var prevSibling = listNode.previousSibling;		
		if(prevSibling == null || prevSibling.tagName == null || prevSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return prevSibling	
	}
	DDList.prototype.isParentBookNode = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		if(this.isNodeItemTypeThis(listNode, ITEMTYPEBOOKOPEN) == true)
			return true;
		else
			return false;
	}
	DDList.prototype.isBookListNode = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		var itemType = this.getNodeItemType(htmlNode);
		if(itemType == ITEMTYPEBOOKOPEN ||
			itemType == ITEMTYPEBOOKCLOSED ||
			itemType == ITEMTYPEPAGE ||
			itemType == ITEMTYPEURL)
			return true;
		else
			return false;
	}
	DDList.prototype.isUrlNode = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		if(this.isNodeItemTypeThis(listNode, ITEMTYPEURL) == true)
			return true;
		else
			return false;	
	}	
	DDList.prototype.setLoadingDisplayInfo = function(iconClass, iconHtml, textClass, textString)
	{
		this.loadingIconClass = iconClass;
		this.loadingIconHtml = iconHtml;
		this.loadingTextClass = textClass;
		this.loadingText = textString;
	}
	DDList.prototype.clearList = function()
	{
		var node = this.rootHtmlNode;
		var childs = node.childNodes;
		var len = childs.length;
		var htmlNode = null;
		var itemType = -1;
		var i=0;
		for(i=len-1; i>=0; i--)
		{
			htmlNode = this.getHtmlNodeFromListNode(childs[i]);
			if(htmlNode != null)
			{
				itemType = this.getNodeItemType(htmlNode);
				if(itemType == ITEMTYPEBOOKOPEN ||
					itemType == ITEMTYPEBOOKCLOSED ||
					itemType == ITEMTYPEPAGE ||
					itemType == ITEMTYPEURL)
					node.removeChild(childs[i]);
			}
		}
	}
	DDList.prototype.insertLoadingMsg = function()
	{
		this.clearList();
		var divLoading = document.createElement('div');
		divLoading.className = LISTITEMCLASS;
		this.setNodeItemType(divLoading, ITEMTYPELOADING);
		var divLoadingIcon = document.createElement('div');
		divLoadingIcon.className = this.loadingIconClass;
		divLoadingIcon.innerHTML = this.loadingIconHtml;
		divLoading.appendChild(divLoadingIcon);
		
		var divLoadingTxt = document.createElement('div');
		divLoadingTxt.className = this.loadingTextClass;
		divLoadingTxt.innerHTML = this.loadingText;
		divLoading.appendChild(divLoadingTxt);
		
		this.rootHtmlNode.appendChild(divLoading);
	}
	DDList.prototype.removeLoadingMsg = function()
	{
		var divLoading = this.getHtmlChildNode(this.rootHtmlNode, "div", ITEMTYPELOADING);
		this.rootHtmlNode.removeChild(divLoading);		
	}
	DDList.prototype.setSelectedListNode = function(listNode)
	{
		if(listNode == null)
			return;
		var nodeClass = "";
		if(this.isParentBookNode(listNode))
			nodeClass = this.openBookClassSelected + " " + UNSELECTABLECLASS;
		else if(this.isBookListNode(listNode))
			nodeClass = this.closedBookClassSelected + " " + UNSELECTABLECLASS;
		else
			nodeClass = this.pageClassSelected + " " + UNSELECTABLECLASS;
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		htmlNode.className = nodeClass;
		
		this.selectedListNode = listNode;
		this.hoveredListNode = listNode;
		
		this.updateSelectedIcon(listNode);
	}	
	DDList.prototype.hoverNode = function(htmlNode, hoverClass)
	{
		if(this.hoveredListNode != null)
		{
			var htmNode = this.getHtmlNodeFromListNode(this.hoveredListNode);
			if(isTouchDevice())
			{
				fireEvent(htmNode, 'touchend');
				fireEvent(htmNode, 'touchmove');
			}
			else
				fireEvent(htmNode, 'mouseout');
		}
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		if(listNode != this.selectedListNode)
			htmlNode.className = hoverClass + " " + UNSELECTABLECLASS;
		this.hoveredListNode = listNode;
		this.updateHoverIcon(listNode);
	}
	DDList.prototype.focusHoveredNode = function()
	{
		if(this.hoveredListNode == null)
			this.hoveredListNode = this.getFirstListNode();
		if(this.hoveredListNode != null)
		{
			var htmlNode = this.getHtmlNodeFromListNode(this.hoveredListNode);
			fireEvent(htmlNode, 'mouseover');
		}
	}
	DDList.prototype.blurHoveredNode = function()
	{
		if(this.hoveredListNode != null)
		{
			var htmlNode = this.getHtmlNodeFromListNode(this.hoveredListNode);
			if(isTouchDevice())
			{
				fireEvent(htmNode, 'touchend');
				fireEvent(htmNode, 'touchmove');
			}
			else
				fireEvent(htmlNode, 'mouseout');
		}
	}
	DDList.prototype.hoverOutNode = function(htmlNode, normalClass)
	{
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		if(listNode != this.selectedListNode)
		{
			htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
			this.updateNormalIcon(listNode);
		}
	}	
	DDList.prototype.clickNode = function(htmlNode, clickClass, normalClass, url)
	{
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		if(listNode == this.selectedListNode)
			return;
		htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
		document.location = url;
	}
	DDList.prototype.clickBook = function(htmlNode)
	{
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		var id = listNode.getAttribute('id');
		var url = this.prepareUrl(id);
		document.location.href = url;
	}
	
	DDList.prototype.loadBook = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		var src = listNode.getAttribute(DATASRC);
		if(src == 'undefined')
		{
			alert(this.errorMsg);
			return;
		}	
		var rootRelPath = listNode.getAttribute(DATAPATH);
		if(rootRelPath == null)
			rootRelPath = this.rootRelPath;
				
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		htmlNode.className = this.openBookClass + " " + UNSELECTABLECLASS;
		this.insertLoadingMsg();
		var parentId = listNode.getAttribute('id');
		this.curParentId = parentId;
		var objContext = new listContext(0,0,0);
			objContext.rootRelPath = rootRelPath;
		var parentBkInfo = new parentBookInfo(rootRelPath, src);
		var srcPath = rootRelPath + "/" + this.dataFolder + "/" + src;
		this.parentBookSrcStack.push(parentBkInfo);
		this.updateUpButton();
		xmlJsReader.loadFile(srcPath, callbackCreateList, objContext);
	}	
	DDList.prototype.isNodeItemTypeThis = function(node, type)
	{
		if(this.getNodeItemType(node) == type)
			return true;
		else
			return false;
	
	}
	DDList.prototype.getIconHtml = function(type)
	{
		if(this.iconHtml[type] != "")
			return this.iconHtml[type];
		var html = "";
		var htmlClass = "";
		var htmlStyle = "";
		var dataItemType = "";
		if(this.iconSrc[type] != "")
		{
			if(this.iconClass[type] != "")
				htmlClass = "class='" + this.iconClass[type] + "' ";
			if(this.iconStyle[type] != "")
				htmlStyle = "class='" + this.iconStyle[type] + "' ";
			dataItemType = DATAITEMTYPE + "='" + ITEMTYPEICON + "' ";
			this.iconHtml[type] = "<img src='" + this.iconSrc[type] + "' " + htmlClass + htmlStyle + dataItemType + "/>";
		}
		return this.iconHtml[type];		
	}
	DDList.prototype.updateHoverIcon = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		var type = this.getNodeItemType(htmlNode);
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		if(iconHtmlNode != null)
		{
			if(this.iconHoverSrc[type] != "")
				iconHtmlNode.setAttribute("src", this.iconHoverSrc[type]);	
		}
	}
	DDList.prototype.updateNormalIcon = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		var type = this.getNodeItemType(htmlNode);
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		if(iconHtmlNode != null)
		{
			if(this.iconSrc[type] != "")
				iconHtmlNode.setAttribute("src", this.iconSrc[type]);	
		}
	}
	DDList.prototype.updateSelectedIcon = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		var type = this.getNodeItemType(htmlNode);
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		if(iconHtmlNode != null)
		{
			if(this.iconSelSrc[type] != "")
				iconHtmlNode.setAttribute("src", this.iconSelSrc[type]);	
		}
	}
	DDList.prototype.getNodeItemType = function(node)
	{
		if(node != null && node != 'undefined')
			return node.getAttribute(DATAITEMTYPE);
		else
			return null;
	}
	DDList.prototype.setNodeItemType = function(node, type)
	{
		if(node != null && node != 'undefined')
			node.setAttribute(DATAITEMTYPE, type);	
	}
	DDList.prototype.prepareUrl = function(id)
	{
		var url = "";
		var urlid = '?' + TOCID + '=' + id;
		var bookmark = getUrlBookmark(url);
		if(bookmark != "")
			url = getUrlWithoutParameterAndBookMark(url);
		urlWithId = url + urlid + bookmark;
		return urlWithId;
	}
	
}
function onNodeHover(htmlNode, hoverClass)
{
	gDDList.hoverNode(htmlNode, hoverClass);
}
function onNodeHoverOut(htmlNode, normalClass)
{
	gDDList.hoverOutNode(htmlNode, normalClass);
}
function onNodeClick(htmlNode, clickClass, normalClass, url)
{
	gDDList.clickNode(htmlNode, clickClass, normalClass, url);
}
function onBackBtnHover(backBtnDiv)
{
	backBtnDiv.className = gDDList.upBtnClassHover;
}
function onBackBtnHoverOut(backBtnDiv)
{
	backBtnDiv.className = gDDList.upBtnClass;
}
function onBackBtnClick(backBtnDiv)
{
	backBtnDiv.className = gDDList.upBtnClassClick;
}
function onKeyPress(e)
{
	gDDList.pressKey(e);
}
function onListFocus()
{
	gDDList.focusHoveredNode();
}
function onListBlur()
{
	gDDList.blurHoveredNode();
}
function callbackCreateList(xmlDoc, arg) //Cannot use binding as IE9 does not support it
{
	gDDList.createChildList(xmlDoc, arg);
}
function onBookClick(htmlNode)
{
	gDDList.clickBook(htmlNode);
}
function onClickUp()
{
	gDDList.clickUp();
}

function loadProjData(childRootRelPath, childCommonRootRelPath, objContext)
{
	if(childRootRelPath == "")
		gDDList.loadFromStack(objContext);
	else
	{
		var strChildProjPath = childRootRelPath + "/" + gDDList.dataFolder + "/" + gDDList.rootFile;
		objContext.rootRelPath = childRootRelPath;
		objContext.commonRootRelPath = childCommonRootRelPath;
		xmlJsReader.loadFile(strChildProjPath, callbackCreateList, objContext);
	}
}
function syncToc(prefix, childOrder)
{
	gDDList.sync(prefix, childOrder);
}
function listContext(bkCount, pgCount, childProjOrder)
{
	this.bookCount = bkCount;
	this.pageCount = pgCount;
	this.rootRelPath = null;
	this.childProjOrder = childProjOrder;	
}
function dataInfo(dataXmlNode, index, rootRelPath, commonRootRelPath, bookCount, pageCount)
{
	this.dataXmlNode = dataXmlNode;
	this.index = index;	
	this.rootRelPath = rootRelPath;
	this.commonRootRelPath = commonRootRelPath;
	this.bookCount = bookCount;
	this.pageCount = pageCount;
}
function parentBookInfo(projPath, src)
{
	this.projPath = projPath;
	this.src = src;	
}