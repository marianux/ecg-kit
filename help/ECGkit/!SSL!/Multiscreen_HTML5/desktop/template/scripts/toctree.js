var	gTopicId = "";
gDocumentURL = "";
function Tree(rootRelPath, dataFolder, rootFile, commonRootRelPath)
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
	
	this.bookChildsClass = "";
	
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

	this.rootFile = rootFile;
	this.rootRelPath = rootRelPath;
	this.commonRootRelPath = commonRootRelPath;
	this.dataFolder = dataFolder;
	this.saveNodesState = true;
	this.rootHtmlNode = null;
	
	this.selectedTreeNode = null;
	this.hoveredTreeNode = null;
	this.syncTree = true;
	
	this.urlId = ""
	this.idParts = null;
	this.idPartsNextIndex = 0;
	this.nextNodeIdToBeSynced = null;
	this.isSyncDone = false;
	this.syncDataLoaded = false;
	
	this.loadingIconClass = "loadingicon";
	this.loadingIconHtml = "";
	this.loadingTextClass = "loadingtext";
	this.loadingText = "";
	this.isSyncingRoot = true;
	
	this.loadStack = new MhStack();

	Tree.prototype.init = function()
	{
		this.rootHtmlNode = document.getElementById(this.rootId);
		if(this.rootHtmlNode.attachEvent)
		{
			this.rootHtmlNode.attachEvent("onkeydown", function(){onKeyPress(event);});
			this.rootHtmlNode.attachEvent("onblur", function(){onTreeBlur();});
			this.rootHtmlNode.attachEvent("onfocus", function(){onTreeFocus();});
		}
		else
		{
			this.rootHtmlNode.setAttribute("onkeydown" , "onKeyPress(event)");	
			this.rootHtmlNode.setAttribute("onblur" , "onTreeBlur()");	
			this.rootHtmlNode.setAttribute("onfocus" , "onTreeFocus()");	
		}
	}	
	Tree.prototype.load = function()
	{
		this.insertLoadingMsg(this.rootHtmlNode);
		var objContext = new treeContext(0,0, this.rootHtmlNode, "", 0);
		xmlJsReader.loadFile(this.rootRelPath + "/" + this.dataFolder + "/" + this.rootFile, callbackCreateTree, objContext);
	}
	Tree.prototype.createChildTree = function(xmlDoc, objContext)
	{
		var dataXmlNode = xmlDoc.getElementsByTagName(DATANODE)[0];
		this.insertTreeItems(dataXmlNode, 0, objContext);
	}
	Tree.prototype.insertTreeItems = function(dataXmlNode, index, objContext)
	{
		var bookCount=objContext.bookCount;
		var pageCount=objContext.pageCount;
		var projOrderStr = objContext.projOrderStr;
		var childProjOrder = objContext.childProjOrder;
		var parentHtmlNode = objContext.parentHtmlNode;
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
				this.insertTreeItem(parentHtmlNode, name, url, bookCount, pageCount, projOrderStr, ITEMTYPEBOOKCLOSED, src, origRootRelPath, origCommonRootRelPath);
			}
			else if(xmlNode.tagName == PAGENODE)
			{
				pageCount++;
				objContext.pageCount = pageCount;
				this.insertTreeItem(parentHtmlNode, name, url, bookCount, pageCount, projOrderStr, ITEMTYPEPAGE);
			}
			else if(xmlNode.tagName == URLNODE)
			{
				pageCount++;
				objContext.pageCount = pageCount;
				this.insertTreeItem(parentHtmlNode, name, url, bookCount, pageCount, projOrderStr, ITEMTYPEURL);
			}
			else if(xmlNode.tagName == PROJNODE)
			{
				++childProjOrder;
				var dInfo = new dataInfo(dataXmlNode, i+1, origRootRelPath, origCommonRootRelPath, bookCount, pageCount, projOrderStr, childProjOrder);
				this.loadStack.push(dInfo);
				objContext.bookCount = 0;
				objContext.pageCount = 0;
				objContext.projOrderStr = projOrderStr + TOCCHILDIDPREFIX + childProjOrder;
				objContext.childProjOrder = 0;
				var ref = xmlNode.getAttribute(REF);
				this.insertChildProjTree(parentHtmlNode, ref, objContext);
				break;
			}
		}
		if(len == i)
			this.loadFromStack(objContext);

	}
	Tree.prototype.loadFromStack = function(objContext)
	{
		if(this.loadStack.isEmpty() != true)
		{
			var dInfo =	this.loadStack.pop();
			objContext.rootRelPath = dInfo.rootRelPath;
			objContext.commonRootRelPath = dInfo.commonRootRelPath;
			objContext.bookCount = dInfo.bookCount;
			objContext.pageCount = dInfo.pageCount;
			objContext.projOrderStr = dInfo.projOrderStr;
			objContext.childProjOrder = dInfo.childProjOrder;
			this.insertTreeItems(dInfo.dataXmlNode, dInfo.index, objContext);
		}
		else
		{
			var parentHtmlNode = objContext.parentHtmlNode;
			if(parentHtmlNode != this.rootHtmlNode)
			{
				var treeNode = this.getTreeNodeFromHtmlNode(parentHtmlNode);
				this.toggleBookNode(treeNode);
				this.isSyncingRoot = false;
			}
			this.removeLoadingMsg(parentHtmlNode);
			if(this.syncTree == true)
			{
				if(this.isSyncDone == false)
				{
					if(this.syncDataLoaded == false)
					{
						this.syncDataLoaded = true;
						loadParentDataForSyncing(gCommonRootRelPath, SCR_PARENT_TOCSYNC);
					}
					else
						this.sync("", "");
				}
			}
		}
	}
	Tree.prototype.insertTreeItem = function(parentHtmlNode, name, url, bookCount, pageCount, projOrderStr, itemType, src, rootRelPath, commonRootRelPath)
	{
		var treeNode = document.createElement("div");
		treeNode.setAttribute('class', TREEITEMCLASS);
		var classNormal = "";
		var classHover = "";
		var classClick = "";
		var html = "";
		var dot = "";
		var parentId = "";
		if(parentHtmlNode != this.rootHtmlNode)
		{
			dot = BOOKDELIM;
			var parentTreeNode = this.getTreeNodeFromHtmlNode(parentHtmlNode);
			parentId = parentTreeNode.getAttribute('id');
		}
		else
			parentId = parentHtmlNode.getAttribute('id');
		
		if(ITEMTYPEBOOKCLOSED == itemType)
		{
			if(src != null && src != '')
				treeNode.setAttribute(DATASRC, src);
			if(url != '' && url != null)
				treeNode.setAttribute(DATAURL, url);
			if(rootRelPath != null)
				treeNode.setAttribute(DATAPATH, rootRelPath);
			if(commonRootRelPath != null)
				treeNode.setAttribute(DATAROOTPATH, commonRootRelPath);
			treeNode.setAttribute('id', parentId + dot + bookCount + projOrderStr);
			classNormal = this.closedBookClass;
			if(url == null)
				classNormal = classNormal + " " + UNCLICKABLECLASS;
			classHover = this.closedBookClassHover;
			classClick = this.closedBookClassClick;
			html = this.bookClosedHtml;
		}
		else if(ITEMTYPEPAGE == itemType)
		{
			treeNode.setAttribute('id', parentId + dot + bookCount + PAGEDELIM + pageCount + projOrderStr);
			classNormal = this.pageClass;
			classHover = this.pageClassHover;
			classClick = this.pageClassClick;
			html = this.pageHtml;
		}
		else if(ITEMTYPEURL == itemType)
		{
			classNormal = this.urlClass;
			classHover = this.urlClassHover;
			classClick = this.urlClassClick;	
			html = this.urlHtml;
		}
		parentHtmlNode.appendChild(treeNode);
		this.insertChildHtmlNode(treeNode, name, itemType, html, classNormal, classHover, classClick, url);
	}
	Tree.prototype.insertChildHtmlNode = function(treeNode, name, itemType, html, classNormal, classHover, classClick, url)
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
			treeNode.insertBefore(anchorNode, treeNode.firstChild);
		}
		else
			treeNode.insertBefore(htmlNode, treeNode.firstChild);
		var newUrl = this.addEventsToNode(htmlNode, classNormal, classHover, classClick, url);
		if(bAddAnchor)
		{
			if(newUrl == "")
				url = "#";
			else 
				url = newUrl;
			if(itemType == ITEMTYPEBOOKCLOSED || itemType == ITEMTYPEBOOKOPEN)
			{
				var curPath = _getPath(document.location.href);
				var absUrl = _getFullPath(curPath, url);
				if(absUrl == document.location.href)
					this.addBookEventsToNode(htmlNode, itemType);
				else
					anchorNode.setAttribute("href", url);
			}
			else
				anchorNode.setAttribute("href", url);
		}
		else
			this.addBookEventsToNode(htmlNode, itemType);
	}
	Tree.prototype.insertChildProjTree = function(parentHtmlNode, ref, objContext)
	{
		var commonRootRelPath = objContext.commonRootRelPath;
		if(commonRootRelPath == null)
			commonRootRelPath = this.commonRootRelPath;
		var strChildProjPath = commonRootRelPath + "/" + ref;
		loadScreenData(strChildProjPath, SCR_CHILD_TOC, objContext);	
	}
	Tree.prototype.addEventsToNode = function(htmlNode, classNormal, classHover, classClick, url)
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
		

		if(url == null || url =='')
			return "";
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		var parentId = treeNode.getAttribute('id');
		var newUrl = "";
		if(parentId != null && parentId != '')
		{
			var rootIdLen = this.rootId.length;
			parentId = parentId.substring(rootIdLen);
				
			var urlid = "";
			if(url.indexOf("?") != -1)
				urlid = '&' + TOCID + '=' + parentId;
			else
				urlid = '?' + TOCID + '=' + parentId;
			var bookmark = getUrlBookmark(url);
			if(bookmark != "")
			{
				var urlWithoutBookmark = getUrlWithoutBookmark(url);
				newUrl = urlWithoutBookmark + urlid + bookmark;
			}
			else
				newUrl = url;
		}
		else
			newUrl = url;
		/*
		if(htmlNode.attachEvent)
			htmlNode.attachEvent('onclick', function(){onNodeClick(htmlNode, classClick, classNormal, newUrl)});
		else
			htmlNode.setAttribute("onclick", "onNodeClick(this, '" + classClick + "', '" + classNormal + "', '" + newUrl + "')");
		*/
		return newUrl;
	}
	Tree.prototype.addBookEventsToNode = function(htmlNode, itemType)
	{
		if(htmlNode.attachEvent)
		{
			if(itemType == ITEMTYPEBOOKCLOSED)
				htmlNode.attachEvent('onclick', function(){onNodeExpand(event, htmlNode);});
			else if(itemType == ITEMTYPEBOOKOPEN)
				htmlNode.attachEvent('onclick', function(){onNodeCollapse(event, htmlNode);});
		}
		else
		{
			if(itemType == ITEMTYPEBOOKCLOSED)
				htmlNode.setAttribute("onclick", "onNodeExpand(event, this)");
			else if(itemType == ITEMTYPEBOOKOPEN)
				htmlNode.setAttribute("onclick", "onNodeCollapse(event, this)");
		}	
	}
	Tree.prototype.toggleBookNode = function(treeNode)
	{
		var classNormal = "";
		var classHover = "";
		var classClick = "";
		var bookHtml = "";
		var curBookHtmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		var itemType = -1;
		
		if(this.isBookOpenTreeNode(treeNode))
		{
			if(treeNode == this.selectedTreeNode)
				classNormal = this.closedBookClassSelected;
			else
				classNormal = this.closedBookClass;
			classHover = this.closedBookClassHover;
			classClick = this.closedBookClassClick;
			bookHtml = this.bookClosedHtml;	
			itemType = ITEMTYPEBOOKCLOSED;
		}
		else
		{
			if(treeNode == this.selectedTreeNode)
				classNormal = this.openBookClassSelected;
			else
				classNormal = this.openBookClass;
			classHover = this.openBookClassHover;
			classClick = this.openBookClassClick;
			bookHtml = this.bookOpenHtml;
			itemType = ITEMTYPEBOOKOPEN;
		}

		var name = curBookHtmlNode.getAttribute("title");
		var url = treeNode.getAttribute(DATAURL);
		var hover = false;
		if(treeNode == this.hoveredTreeNode)
			hover = true;
		if(url == null || url == 'undefinied')
			treeNode.removeChild(curBookHtmlNode);
		else
			treeNode.removeChild(curBookHtmlNode.parentNode);
		this.insertChildHtmlNode(treeNode, name, itemType, bookHtml, classNormal, classHover, classClick, url);
		if(treeNode == this.selectedTreeNode)
			this.updateSelectedIcon(treeNode);
		var bookHtmlNode = this.getHtmlNodeFromTreeNode(treeNode);
	}
	Tree.prototype.expandTreeNode = function(treeNode)
	{
		var src = treeNode.getAttribute(DATASRC);
		if(src == 'undefined')
		{
			alert(this.errorMsg);
			return;
		}
		if(this.isBookNodeLoaded(treeNode))
		{
			this.toggleBookNode(treeNode);
			var bookChildsNode = treeNode.lastChild;
			bookChildsNode.style.display = "block";
		}
		else
		{
			var rootRelPath = treeNode.getAttribute(DATAPATH);
			if(rootRelPath == null)
				rootRelPath = this.rootRelPath;
			var commonRootRelPath = treeNode.getAttribute(DATAROOTPATH);
			if(commonRootRelPath == null)
				commonRootRelPath = this.commonRootRelPath;
			var bookChildsNode = document.createElement('div');
			bookChildsNode.className = this.bookChildsClass;
			bookChildsNode.style.display = "block";
			this.setNodeItemType(bookChildsNode, ITEMTYPEBOOKCHILDS);
			treeNode.appendChild(bookChildsNode);
			this.insertLoadingMsg(bookChildsNode);
			var id = treeNode.getAttribute("id");
			var projOrderStr = this.getProjOrderStr(id);
			var objContext = new treeContext(0,0, bookChildsNode, projOrderStr, 0);
			objContext.rootRelPath = rootRelPath;
			objContext.commonRootRelPath = commonRootRelPath;
			xmlJsReader.loadFile(rootRelPath + "/" + this.dataFolder + "/" + src, callbackCreateTree, objContext);
		}
	}
	Tree.prototype.collapseTreeNode = function(treeNode)
	{
		this.toggleBookNode(treeNode);
		var bookChildsNode = treeNode.lastChild;
		bookChildsNode.style.display = "none";
	}
	Tree.prototype.sync = function(childPrefix, childOrder)
	{
		if(this.isSyncDone == true)
			return;
		if(childPrefix != "")
			childPrefix += ".";
		if(this.isSyncingRoot)
		{
			this.isSyncingRoot = false;
			this.urlId = getUrlParameter(TOCID, gDocumentURL);
			if(this.urlId == '')
			{
				if(gTopicId != null && gTopicId != 'undefined' && gTopicId != '')
				{
					var dataId = this.rootId;
					gTopicId = gTopicId.split(BOOKDELIM).join(childOrder + BOOKDELIM);
					this.urlId = dataId + childPrefix + gTopicId + childOrder;
				}
				else
				{
					this.isSyncDone = true;
					return;
				}
			}
			else
				this.urlId = this.rootId + this.urlId;
			this.idParts = this.urlId.split(BOOKDELIM);
			this.nextNodeIdToBeSynced = this.idParts[this.idPartsNextIndex];
		}
		var treeNode = this.getTreeNodeById(this.nextNodeIdToBeSynced);
		++this.idPartsNextIndex;
		if(this.idPartsNextIndex < this.idParts.length)
		{
			this.expandTreeNode(treeNode);
			this.nextNodeIdToBeSynced = this.nextNodeIdToBeSynced + BOOKDELIM + this.idParts[this.idPartsNextIndex];
		}
		else
		{
			this.isSyncDone = true;
			this.setSelectedTreeNode(treeNode);
			if(this.isBookClosedTreeNode(treeNode))
				this.expandTreeNode(treeNode);
		}
	}
	Tree.prototype.pressKey = function(e)
	{
		if (e.preventDefault)
            e.preventDefault();
		var kCode = 0;
		if(e.keyCode)
			kCode = e.keyCode;
		else
			kCode = e.which;
		var treeNode = null;
		var htmlNode = null;
		var event = "";
		if(kCode == 38)
		{
			treeNode = this.getPreviousTreeItem(this.hoveredTreeNode);
			htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
			event = "mouseover";
		}
		else if(kCode == 40)
		{
			treeNode = this.getNextTreeItem(this.hoveredTreeNode);
			htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
			event = "mouseover";
		}
		else if(kCode == 13 || kCode == 32)
		{
			treeNode = this.hoveredTreeNode;
			htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
			event = "click";
		}
		else if(kCode == 39)
		{
			if(this.isBookClosedTreeNode(this.hoveredTreeNode))
			{
				treeNode = this.hoveredTreeNode;
				htmlNode = this.getIconHtmlNodeFromTreeNode(treeNode);
				event = "click";
			}
		}
		else if(kCode == 37)
		{
			if(this.isBookOpenTreeNode(this.hoveredTreeNode))
			{
				treeNode = this.hoveredTreeNode;
				htmlNode = this.getIconHtmlNodeFromTreeNode(treeNode);
				event = "click";
			}
			else
			{
				treeNode = this.getParentTreeNode(this.hoveredTreeNode)
				htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
				event = "mouseover";
			}
		}
		if(htmlNode != null)
			fireEvent(htmlNode, event);
	}
	Tree.prototype.getProjOrderStr = function(id)
	{
		var pos = id.lastIndexOf(BOOKDELIM);
		var orderStr = "";
		if(pos != -1)
			orderStr = id.substring(pos);
		else
			orderStr = id;
		pos = orderStr.indexOf(TOCCHILDIDPREFIX);
		if(pos != -1)
			orderStr = orderStr.substring(pos);
		else
			orderStr = "";
		return orderStr;
	}
	Tree.prototype.getTreeNodeById = function(id)
	{
		return document.getElementById(id);
	}
	Tree.prototype.getHtmlNodeFromTreeNode = function(treeNode)
	{
		if(treeNode == null)
			return null;
		var childHtmlNodes = null;
		var childs = treeNode.childNodes;
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
			childHtmlNodes = treeNode.getElementsByTagName("div");
		var htmlNode = null;
		htmlNode = childHtmlNodes[0];
		return htmlNode;
	}
	Tree.prototype.getTreeNodeFromHtmlNode = function(htmlNode)
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
	Tree.prototype.getIconHtmlNodeFromTreeNode = function(treeNode)
	{
		if(treeNode == null)
			return null;
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode)
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		return iconHtmlNode;
	}	
	Tree.prototype.getTreeNodeFromIconHtmlNode = function(iconHtmlNode)
	{
		if(iconHtmlNode == null)
			return null;
		var parentHtmlNode = iconHtmlNode.parentNode;
		while(parentHtmlNode != null && parentHtmlNode != 'undefined' && parentHtmlNode != this.rootHtmlNode)
		{
			var itemType = this.getNodeItemType(parentHtmlNode);
			if(itemType == ITEMTYPEBOOKOPEN || itemType == ITEMTYPEBOOKCLOSED || itemType == ITEMTYPEPAGE || itemType == ITEMTYPEURL)
				break;
			parentHtmlNode = parentHtmlNode.parentNode;
		}
		if(parentHtmlNode != this.rootHtmlNode)
			return this.getTreeNodeFromHtmlNode(parentHtmlNode);
		else
			return null;
	}
	Tree.prototype.getHtmlChildNode = function(node, tag, type)
	{
		if(tag == "" || tag == 'undefined' || node == null)
			return null;
		var nodeChilds = node.getElementsByTagName(tag);
		var len = nodeChilds.length;
		var i=0;
		for(i=0; i<len; i++)
		{
			if(this.isNodeItemTypeThis(nodeChilds[i],type))
			 return nodeChilds[i];		
		}
		return null;
	}
	Tree.prototype.getFirstTreeNode = function()
	{
		var treeNodes = this.rootHtmlNode.getElementsByTagName("div");
		if(treeNodes.length > 0)
			return treeNodes[0];
	}
	Tree.prototype.getNextTreeItem = function(treeNode)
	{
		var nextNode = null;
		if(treeNode == null || treeNode == 'undefined')
			return null;
		if(this.isBookOpenTreeNode(treeNode))
			nextNode = this.getFirstChildNode(treeNode);
		if(nextNode == null)
			nextNode = treeNode.nextSibling; 
		if(nextNode == null)
		{
			var parentBookNode = treeNode;

			while(nextNode == null)
			{
				parentBookNode = this.getParentTreeNode(parentBookNode);
				if(parentBookNode == null)
					break;
				nextNode = this.getNextSiblingNode(parentBookNode);
			}
		}
		return nextNode;	
	}
	Tree.prototype.getPreviousTreeItem = function(treeNode)
	{
		var prevNode = this.getPreviousSiblingNode(treeNode);
		if(prevNode != null)
		{
			var lastChildNode = prevNode;
			while(lastChildNode != null)
			{
				prevNode = lastChildNode;
				lastChildNode = this.getLastChildNode(lastChildNode);
			}			
		}
		if(prevNode == null)
			prevNode = this.getParentTreeNode(treeNode);
		return prevNode;	
	}
	Tree.prototype.getFirstChildNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var bookChildsNode = null;
		if(this.isBookOpenTreeNode(treeNode))
			bookChildsNode = treeNode.getElementsByTagName("div")[1];
		if(bookChildsNode != null)
			return bookChildsNode.firstChild;
		else
			return null;
	}
	Tree.prototype.getLastChildNode = function(treeNode)
	{
		var lastChild = null;
		if(treeNode == null || treeNode == 'undefined')
			return null;
		lastChild = this.getFirstChildNode(treeNode);
		nextSibling = lastChild;
		while(nextSibling != null)
		{
			lastChild = nextSibling;
			nextSibling = this.getNextSiblingNode(nextSibling);
		}
		return lastChild;
	}
	Tree.prototype.getParentTreeNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var bookChildsNode = treeNode.parentNode;
		if(bookChildsNode != null && this.isNodeItemTypeThis(bookChildsNode,ITEMTYPEBOOKCHILDS) == false)
			return null;
		else
			return bookChildsNode.parentNode;
	}
	Tree.prototype.getNextSiblingNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var nextSibling = treeNode.nextSibling;
		if(nextSibling == null || nextSibling.tagName == null || nextSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return nextSibling
	}
	Tree.prototype.getPreviousSiblingNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var prevSibling = treeNode.previousSibling;		
		if(prevSibling == null || prevSibling.tagName == null || prevSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return prevSibling
	}
	Tree.prototype.isBookOpenTreeNode = function(treeNode)
	{
		if(treeNode == null)
			return false;
		var src = treeNode.getAttribute(DATASRC);
		if(src == null || src == '')
			return false;
		else
		{
			var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
			if(this.isNodeItemTypeThis(htmlNode,ITEMTYPEBOOKOPEN))
				return true;
			else
				return false;
		}
	}
	Tree.prototype.isBookClosedTreeNode = function(treeNode)
	{
		if(treeNode == null)
			return false;
		var src = treeNode.getAttribute(DATASRC);
		if(src == null || src == '')
			return false;
		else
		{
			var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
			if(this.isNodeItemTypeThis(htmlNode,ITEMTYPEBOOKCLOSED))
				return true;
			else
				return false;
		}
	}	
	Tree.prototype.isUrlNode = function(treeNode)
	{
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		if(this.isNodeItemTypeThis(htmlNode, ITEMTYPEURL))
			return true;
		else
			return false;	
	}
	Tree.prototype.setLoadingDisplayInfo = function(iconClass, iconHtml, textClass, textString)
	{
		this.loadingIconClass = iconClass;
		this.loadingIconHtml = iconHtml;
		this.loadingTextClass = textClass;
		this.loadingText = textString;
	}
	Tree.prototype.insertLoadingMsg = function(htmlNode)
	{
		var divLoading = document.createElement('div');
		divLoading.className = TREEITEMCLASS;
		this.setNodeItemType(divLoading, ITEMTYPELOADING);
		var divLoadingImg = document.createElement('div');
		divLoadingImg.className = this.loadingIconClass;
		divLoadingImg.innerHTML = this.loadingIconHtml ; 
		divLoading.appendChild(divLoadingImg);
		var divLoadingTxt = document.createElement('div');
		divLoadingTxt.className = this.loadingTextClass;
		divLoadingTxt.innerHTML = this.loadingText;
		divLoading.appendChild(divLoadingTxt);
		htmlNode.appendChild(divLoading);
	}
	Tree.prototype.removeLoadingMsg = function(htmlNode)
	{
		var divLoading = this.getHtmlChildNode(htmlNode, "div", ITEMTYPELOADING);
		htmlNode.removeChild(divLoading);	
	}
	Tree.prototype.setSelectedTreeNode = function(treeNode)
	{
		if(treeNode == null)
			return;
		var nodeClass = "";
		if(this.isBookClosedTreeNode(treeNode))
			nodeClass = this.closedBookClassSelected + " " + UNSELECTABLECLASS;
		else if(this.isBookOpenTreeNode(treeNode))
			nodeClass = this.openBookClassSelected + " " + UNSELECTABLECLASS;	
		else
			nodeClass = this.pageClassSelected + " " + UNSELECTABLECLASS;
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		htmlNode.className = nodeClass;
		
		this.selectedTreeNode = treeNode;
		this.hoveredTreeNode = treeNode;
		
		this.updateSelectedIcon(treeNode);
	}
	Tree.prototype.hoverNode = function(htmlNode, hoverClass)
	{
		if(this.hoveredTreeNode != null)
		{
			var htmNode = this.getHtmlNodeFromTreeNode(this.hoveredTreeNode);
			if(isTouchDevice())
			{
				fireEvent(htmNode, 'touchend');
				fireEvent(htmNode, 'touchmove');
			}
			else
				fireEvent(htmlNode, 'mouseout');
		}
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		if(treeNode != this.selectedTreeNode)
			htmlNode.className = hoverClass + " " + UNSELECTABLECLASS;
		this.hoveredTreeNode = treeNode;
		if(treeNode != this.selectedTreeNode)
			this.updateHoverIcon(treeNode);
	}
	Tree.prototype.focusHoveredNode = function()
	{
		if(this.hoveredTreeNode == null)
			this.hoveredTreeNode = this.getFirstTreeNode();
		if(this.hoveredTreeNode != null)
		{
			var htmlNode = this.getHtmlNodeFromTreeNode(this.hoveredTreeNode);
			fireEvent(htmlNode, 'mouseover');
		}
	}
	Tree.prototype.blurHoveredNode = function()
	{
		if(this.hoveredTreeNode != null)
		{
			var htmlNode = this.getHtmlNodeFromTreeNode(this.hoveredTreeNode);
			if(isTouchDevice())
			{
				fireEvent(htmNode, 'touchend');
				fireEvent(htmNode, 'touchmove');
			}
			else
				fireEvent(htmlNode, 'mouseout');
		}
	}
	Tree.prototype.hoverOutNode = function(htmlNode, normalClass)
	{
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		if(treeNode != this.selectedTreeNode)
		{
			htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
			this.updateNormalIcon(treeNode);
		}
	}	
	Tree.prototype.clickNode = function(htmlNode, clickClass, normalClass, url)
	{
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		if(treeNode == this.selectedTreeNode)
			return;
		htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
		document.location = url;
	}
	Tree.prototype.isBookNodeLoaded = function(treeNode)
	{
		if(treeNode.getElementsByTagName("div").length > 2)
			return true;
		else
			return false;	
	}
	Tree.prototype.isNodeItemTypeThis = function(node, type)
	{
		if(this.getNodeItemType(node) == type)
			return true;
		else
			return false;
	
	}
	Tree.prototype.getIconHtml = function(type)
	{
		if(this.iconHtml[type] != "")
			return this.iconHtml[type];
		var html = "";
		var htmlClass = "";
		var htmlStyle = "";
		var htmlClickEvt = "";
		var dataItemType = "";
		if(this.iconSrc[type] != "")
		{
			if(this.iconClass[type] != "")
				htmlClass = "class='" + this.iconClass[type] + "' ";
			if(this.iconStyle[type] != "")
				htmlStyle = "class='" + this.iconStyle[type] + "' ";
			if(ITEMTYPEBOOKCLOSED == type)
				htmlClickEvt = "onclick='onNodeExpand(event, this)' ";
			else if(ITEMTYPEBOOKOPEN == type)
				htmlClickEvt = "onclick='onNodeCollapse(event, this)' ";
			dataItemType = DATAITEMTYPE + "='" + ITEMTYPEICON + "' ";
			this.iconHtml[type] = "<img src='" + this.iconSrc[type] + "' " + htmlClass + htmlStyle + htmlClickEvt + dataItemType + "/>";
		}
		return this.iconHtml[type];		
	}	
	Tree.prototype.updateHoverIcon = function(treeNode)
	{
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		var type = this.getNodeItemType(htmlNode);
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		if(iconHtmlNode != null)
		{
			if(this.iconHoverSrc[type] != "")
				iconHtmlNode.setAttribute("src", this.iconHoverSrc[type]);	
		}
	}
	Tree.prototype.updateNormalIcon = function(treeNode)
	{
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		var type = this.getNodeItemType(htmlNode);
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		if(iconHtmlNode != null)
		{
			if(this.iconSrc[type] != "")
				iconHtmlNode.setAttribute("src", this.iconSrc[type]);	
		}
	}
	Tree.prototype.updateSelectedIcon = function(treeNode)
	{
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		var type = this.getNodeItemType(htmlNode);
		var iconHtmlNode = this.getHtmlChildNode(htmlNode, "img", ITEMTYPEICON);
		if(iconHtmlNode != null)
		{
			if(this.iconSelSrc[type] != "")
				iconHtmlNode.setAttribute("src", this.iconSelSrc[type]);	
		}
	}
	Tree.prototype.getNodeItemType = function(node)
	{
		if(node != null && node != 'undefined')
			return node.getAttribute(DATAITEMTYPE);
		else
			return null;
	}
	Tree.prototype.setNodeItemType = function(node, type)
	{
		if(node != null && node != 'undefined')
			node.setAttribute(DATAITEMTYPE, type);	
	}
}

function onNodeHover(htmlNode, hoverClass)
{
	gTree.hoverNode(htmlNode, hoverClass);
}
function onNodeHoverOut(htmlNode, normalClass)
{
	gTree.hoverOutNode(htmlNode, normalClass);
}
function onNodeClick(htmlNode, clickClass, normalClass, url)
{
	gTree.clickNode(htmlNode, clickClass, normalClass, url);
}
function onKeyPress(e)
{
	gTree.pressKey(e);
}
function onTreeFocus()
{
	gTree.focusHoveredNode();
}
function onTreeBlur()
{
	gTree.blurHoveredNode();
}
function callbackCreateTree(xmlDoc, arg) //Cannot use binding as IE9 does not support it
{
	gTree.createChildTree(xmlDoc, arg);
}
function onNodeExpand(e, htmlNode)
{
	var evt = e || window.event; // IE compatibility
	if(evt.preventDefault)
	{
		evt.preventDefault();  
	}
	else
	{
		evt.returnValue = false;  
		evt.cancelBubble=true;  
	}
	var treeNode = null;
	if(gTree.isNodeItemTypeThis(htmlNode, ITEMTYPEBOOKCLOSED))
		treeNode = gTree.getTreeNodeFromHtmlNode(htmlNode);
	else if(gTree.isNodeItemTypeThis(htmlNode, ITEMTYPEICON))
		treeNode = gTree.getTreeNodeFromIconHtmlNode(htmlNode);
	
	if(treeNode)
		gTree.expandTreeNode(treeNode);
}
function onNodeCollapse(e, htmlNode)
{
	var evt = e || window.event; // IE compatibility
	if(evt.preventDefault)
	{
		evt.preventDefault();  
	}
	else
	{
		evt.returnValue = false;  
		evt.cancelBubble=true;  
	}
	var treeNode = null;
	if(gTree.isNodeItemTypeThis(htmlNode, ITEMTYPEBOOKOPEN))
		treeNode = gTree.getTreeNodeFromHtmlNode(htmlNode);
	else if(gTree.isNodeItemTypeThis(htmlNode, ITEMTYPEICON))
		treeNode = gTree.getTreeNodeFromIconHtmlNode(htmlNode);

	if(treeNode)
		gTree.collapseTreeNode(treeNode);
}
function loadProjData(childRootRelPath, childCommonRootRelPath, objContext)
{
	if(childRootRelPath == "")
		gTree.loadFromStack(objContext);
	else
	{
		var strChildProjPath = childRootRelPath + "/" + gTree.dataFolder + "/" + gTree.rootFile;
		objContext.rootRelPath = childRootRelPath;
		objContext.commonRootRelPath = childCommonRootRelPath;
		xmlJsReader.loadFile(strChildProjPath, callbackCreateTree, objContext);	
	}
}
function syncToc(prefix, childOrder)
{
	gTree.sync(prefix, childOrder);
}
function treeContext(bkCount, pgCount, parentHtmlNode, projOrderStr, childProjOrder)
{
	this.bookCount = bkCount;
	this.pageCount = pgCount;
	this.parentHtmlNode = parentHtmlNode;
	this.rootRelPath = null;
	this.commonRootRelPath = null;
	this.projOrderStr = projOrderStr;
	this.childProjOrder = childProjOrder;
}
function dataInfo(dataXmlNode, index, rootRelPath, commonRootRelPath, bookCount, pageCount, projOrderStr, childProjOrder)
{
	this.dataXmlNode = dataXmlNode;
	this.index = index;	
	this.rootRelPath = rootRelPath;
	this.commonRootRelPath = commonRootRelPath;
	this.bookCount = bookCount;
	this.pageCount = pageCount;
	this.projOrderStr = projOrderStr;
	this.childProjOrder = childProjOrder;
}