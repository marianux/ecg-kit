function GloRootFileXmlObject(xmlDoc, i)
{
	this.xmlDoc = xmlDoc;
	this.nodeIndex = i;
}
function GloChunkXmlObject(xmlDoc, i, len)
{
	this.xmlDoc = xmlDoc;
	this.nodeIndex = i;
	this.length = len;
}
function TermInfoObject(node, parentObj)
{
	this.node = node;
	this.parentObj = parentObj;
}
function GloList(gloRootPathsArr, dataFolder, rootFile)
{
	this.rootPathsArr = gloRootPathsArr;
	this.rootFilesXmlArr = new Array;
	this.curChunksXmlArr = new Array;
	this.nextChunkIndex = 0;
	this.curTermName = "";
	this.curChunksToBeMerged = new Array;
	this.curCategory = "";

	this.errorMsg = "";
	this.termClass = "";
	this.termStyle = "";
	this.termClassHover = "";
	this.termClassClick = "";
	this.defClass = "";
	this.defStyle = "";
	this.defClassHover = "";
	this.defClassClick = "";
	this.categoryClass = "";
	this.categoryStyle = "";
	this.termHtml = "";
	this.defHtml = "";
	this.categoryHtml = "";
	
	this.filterBoxId = "";
	
	this.dataFolder = dataFolder;
	this.saveNodesState = true;
	this.rootHtmlNode = null;
	
	this.rootFile = rootFile;
	this.selectedListNode = null;
	this.hoveredListNode = null;
	this.syncTree = true;
	
	this.urlId = "";
	this.idParts = null;
	this.idPartsNextIndex = 0;
	this.nextNodeIdToBeSynced = null;
	this.doSyncNeeded = true;
	
	this.loadingIconClass = "loadingicon";
	this.loadingIconHtml = "";
	this.loadingTextClass = "loadingtext";
	this.loadingText = "";
	
	GloList.prototype.init = function()
	{
		this.rootHtmlNode = document.getElementById(this.rootId);
		if(this.rootHtmlNode.attachEvent)
		{
			this.rootHtmlNode.attachEvent("onkeydown", function(){onGloListKeyPress(event);});
			this.rootHtmlNode.attachEvent("onblur", function(){onGloListBlur();});
			this.rootHtmlNode.attachEvent("onfocus", function(){onGloListFocus();});
		}
		else
		{
			this.rootHtmlNode.setAttribute("onkeydown" , "onGloListKeyPress(event)");	
			this.rootHtmlNode.setAttribute("onblur" , "onGloListBlur()");	
			this.rootHtmlNode.setAttribute("onfocus" , "onGloListFocus()");	
		}
		if(gbIE5)
			readSetting(RHCSHMODE, callbackGloCSHModeRead);
	}
	GloList.prototype.load = function()
	{
		this.insertLoadingMsg(this.rootHtmlNode);
		var i=0;
		xmlJsReader.loadFile(this.rootPathsArr[0] + "/" + this.dataFolder + "/" + this.rootFile, callbackGloRootFileLoaded, i);
	}
	GloList.prototype.loadRootFiles = function(xmlDoc, i)
	{
		var rootFileXmlObj = new GloRootFileXmlObject(xmlDoc, 0);
		this.rootFilesXmlArr[i] = rootFileXmlObj;
		chunkXmlObj = new GloChunkXmlObject(null, 0, 0);
		this.curChunksXmlArr[i] = chunkXmlObj;
		++i;
		var len = this.rootPathsArr.length;
		if(i < len)
			xmlJsReader.loadFile(this.rootPathsArr[i] + "/" + this.dataFolder + "/" + this.rootFile, callbackGloRootFileLoaded, i);
		else
			this.mergeTerms();
	}
	GloList.prototype.mergeTerms = function()
	{
		var chunkXmlObj = this.curChunksXmlArr[this.nextChunkIndex];
		if(chunkXmlObj.xmlDoc == null || chunkXmlObj.nodeIndex < 0 || chunkXmlObj.nodeIndex >= chunkXmlObj.length)
		{
			var rootFileXmlObj = this.rootFilesXmlArr[this.nextChunkIndex];
			if(rootFileXmlObj.nodeIndex == -1)
			{
				this.incrementNextChunkIndexAndDoAction();
				return;
			}
			var indexXmlNode = rootFileXmlObj.xmlDoc.getElementsByTagName(GLOSSARYNODE)[0];
			var childNodes = indexXmlNode.getElementsByTagName(CHUNKINFONODE);
			var len = childNodes.length;
			if(rootFileXmlObj.nodeIndex < len)
			{
				var chunkInfoNode = childNodes[rootFileXmlObj.nodeIndex];
				rootFileXmlObj.nodeIndex++;
				var url = chunkInfoNode.getAttribute(URL);
				var chunkPath = this.rootPathsArr[this.nextChunkIndex] + "/" + this.dataFolder + "/" + url;
				xmlJsReader.loadFile(chunkPath, callbackGloChunkLoaded);
			}
			else
			{
				rootFileXmlObj.nodeIndex = -1;
				this.incrementNextChunkIndexAndDoAction();
				return;
			}
			
		}
		else
		{
			var dataNode = chunkXmlObj.xmlDoc.getElementsByTagName(DATANODE)[0];
			var childNodes = getChildElementsByTagName(dataNode,ENTRYNODE);
			var entryNode = childNodes[chunkXmlObj.nodeIndex];
			var name = entryNode.getAttribute(NAME);
			if(this.curTermName == "" || compare(name,this.curTermName) == -1)
			{
				this.curTermName = name;
				this.curChunksToBeMerged.splice(0, this.curChunksToBeMerged.length);
				var termInfoObj = new TermInfoObject(entryNode, chunkXmlObj);
				this.curChunksToBeMerged[0] = termInfoObj;
			}
			else if(name == this.curTermName)
			{
				var termInfoObj = new TermInfoObject(entryNode, chunkXmlObj);
				this.curChunksToBeMerged[this.curChunksToBeMerged.length] = termInfoObj;
			}
			this.incrementNextChunkIndexAndDoAction();
			return;
		}
	}
	GloList.prototype.incrementNextChunkIndexAndDoAction = function()
	{
		this.nextChunkIndex++;
		if(this.nextChunkIndex >= this.rootPathsArr.length)
		{
			this.nextChunkIndex = 0;
			if(this.curChunksToBeMerged.length != 0)
			{
				this.insertTerm(this.rootHtmlNode, this.curChunksToBeMerged);
				this.curChunksToBeMerged.splice(0, this.curChunksToBeMerged.length);
				this.curTermName = "";
			}
			else
			{
				this.removeLoadingMsg(this.rootHtmlNode);
				this.curTermName = "";
				this.filterKeywords(true);
				return;
			}
		}
		this.mergeTerms();
	}
	GloList.prototype.readChunk = function(xmlDoc, arg)
	{
		var chunkXmlObj = this.curChunksXmlArr[this.nextChunkIndex];
		var dataNode = xmlDoc.getElementsByTagName(DATANODE)[0];
		chunkXmlObj.xmlDoc = xmlDoc;
		chunkXmlObj.nodeIndex = 0;
		chunkXmlObj.length = getChildElementsByTagName(dataNode,ENTRYNODE).length;
		this.mergeTerms();
	}
	GloList.prototype.insertTerm = function(parentHtmlNode, chunksArr)
	{
		var classNormal = this.termClass;
		var classHover = this.termClassHover;
		var classClick = this.termClassClick;
		var inlinestyle = this.termStyle;
		
		var defClassNormal = this.defClass;
		var defClassHover = this.defClassHover;
		var defClassClick = this.defClassClick;
		var defInlineStyle = this.defStyle;
		var html = this.termHtml;
		var defHtml = this.defHtml;
		

		var termName = this.curTermName;
		var ch = termName.substring(0, 1).toLocaleUpperCase();
		if(compare(ch,this.curCategory)!=0)
		{
			this.curCategory = ch;
			this.insertCategory(ch);
		}
		
		var listNode = document.createElement("div");
		listNode.setAttribute('class', TREEITEMCLASS);
		if(parentHtmlNode == this.rootHtmlNode)
		{
			var divLoading = this.getLoadingHtmlNode();
			parentHtmlNode.insertBefore(listNode, divLoading);
		}
		else
			parentHtmlNode.appendChild(listNode);
		

		this.insertChildHtmlNode(listNode, termName, ITEMTYPETERM, html, classNormal, classHover, classClick, inlinestyle);
		
		var len = chunksArr.length;
		var def = "";
		var br = "";
		for(var i=0; i<len; i++)
		{
			var termInfoObj = chunksArr[i];
			if(termInfoObj.parentObj != null)
				termInfoObj.parentObj.nodeIndex++;
			var entryNode = termInfoObj.node;
			if(i>0)
				br = "<br />";
			def = def + br + entryNode.getAttribute(VALUE);
		}
		this.insertChildHtmlNode(listNode, def, ITEMTYPEDEF, defHtml, defClassNormal, defClassHover, defClassClick, defInlineStyle);

	}
	GloList.prototype.insertChildHtmlNode = function(parentHtmlNode, name, itemType, html, classNormal, classHover, classClick, style)
	{
		html = html.replace(LINK_NAME_MACRO, name);

		var htmlNode = document.createElement("div");
		htmlNode.className = classNormal + " " + UNSELECTABLECLASS;
		if(style != "")
			htmlNode.style.cssText = style;
		htmlNode.innerHTML = html;

		parentHtmlNode.appendChild(htmlNode);
		this.addEventsToNode(htmlNode, classNormal, classHover, classClick);

		this.setNodeItemType(htmlNode, itemType);
		if(itemType == ITEMTYPETERM)
		{
			htmlNode.setAttribute("title", name);
			this.setNodeTerm(htmlNode, name);
		}
		else if(itemType == ITEMTYPEDEF)
			htmlNode.style.display = "none";
	}
	GloList.prototype.insertCategory = function(ch)
	{
		var listNode = document.createElement("div");
		listNode.setAttribute('class', TREEITEMCLASS);
		
		var categoryElem = document.createElement("div");
		categoryElem.className = this.categoryClass;
		categoryElem.style.cssText = this.categoryStyle;
		categoryElem.innerHTML = this.categoryHtml.replace(LINK_NAME_MACRO, ch);
		this.setNodeItemType(categoryElem, ITEMTYPECATEGORY);
		listNode.appendChild(categoryElem);
		var divLoading = this.getLoadingHtmlNode();
		this.rootHtmlNode.insertBefore(listNode, divLoading);
	}
	GloList.prototype.addEventsToNode = function(htmlNode, classNormal, classHover, classClick, url)
	{
		if(htmlNode.attachEvent)
		{
			if(isTouchDevice())
			{
				htmlNode.attachEvent('ontouchstart', function(){onGloNodeHover(htmlNode, classHover);});
				htmlNode.attachEvent('ontouchend', function(){onGloNodeHoverOut(htmlNode, classHover);});
				htmlNode.attachEvent('ontouchmove', function(){onGloNodeHoverOut(htmlNode, classHover);});
			}
			else
			{
				htmlNode.attachEvent('onmouseout', function(){onGloNodeHoverOut(htmlNode, classNormal);});
				htmlNode.attachEvent('onmouseover', function(){onGloNodeHover(htmlNode, classHover);});
			}
			htmlNode.attachEvent('onclick', function(){onGloNodeClick(htmlNode);});
		}
		else
		{
			if(isTouchDevice())
			{
				htmlNode.setAttribute("ontouchstart", "onGloNodeHover(this,'" + classHover + "')");
				htmlNode.setAttribute("ontouchend", "onGloNodeHoverOut(this,'" + classNormal + "')");
			}
			else
			{
				htmlNode.setAttribute("onmouseout", "onGloNodeHoverOut(this,'" + classNormal + "')");
				htmlNode.setAttribute("onmouseover", "onGloNodeHover(this,'" + classHover + "')");
			}
			htmlNode.setAttribute("onclick", "onGloNodeClick(this)");
		}
	}
	GloList.prototype.filterKeywords = function(bIsOnLoad)
	{
		var filterBox = document.getElementById(this.filterBoxId);
		if(filterBox == null)
			return true;
		var strFilter = filterBox.value;
		var placeholderText = filterBox.getAttribute(DATAPH);
		if(gbIE55 && !gbIE10 && placeholderText == strFilter)
			strFilter = "";	
		strFilter = strFilter.toLocaleLowerCase();

		if(bIsOnLoad && strFilter == "")
			return true;
		var listNodes = this.rootHtmlNode.childNodes;
		var foundDisplayNode = false;
		var itemDisplayed = false;
		var categoryNode = null;
		for (var i = 0; i < listNodes.length; i++)
		{
			listNode = listNodes[i];
			if(listNode.nodeType != JS_TAGTOKEN)
				continue;
			var htmlNode = this.getHtmlNodeFromListNode(listNode);
			var itemType = this.getNodeItemType(htmlNode);
			if(itemType == ITEMTYPETERM)
			{
				var strItemValue = this.getNodeTerm(htmlNode).toLocaleLowerCase();
				if(strItemValue.indexOf(strFilter)>=0)
				{
					listNode.style.display = "block";
					foundDisplayNode = true;
					itemDisplayed = true;
				}
				else
					listNode.style.display = "none";
			}
			else if(itemType == ITEMTYPECATEGORY) 
			{
				if(i>0 && categoryNode != null)
				{
					if(foundDisplayNode == false)
						categoryNode.style.display = "none";
					else
						categoryNode.style.display = "block";
				}
				categoryNode = listNode;
				foundDisplayNode = false;
			}

		}
		if(categoryNode != null)
		{
			if(foundDisplayNode == false)
				categoryNode.style.display = "none";
			else
				categoryNode.style.display = "block";
		}
		return true;
	}
	GloList.prototype.pressKey = function(e)
	{

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
			listNode = this.getPreviousTreeItem(this.hoveredListNode);
			htmlNode = this.getHtmlNodeFromListNode(listNode);
			event = "mouseover";
		}
		else if(kCode == 40)
		{
			listNode = this.getNextTreeItem(this.hoveredListNode);
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
			if(this.isBookClosedListNode(this.hoveredListNode))
			{
				listNode = this.hoveredListNode;
				htmlNode = this.getIconHtmlNodeFromListNode(listNode);
				event = "click";
			}
		}
		else if(kCode == 37)
		{
			if(this.isBookOpenListNode(this.hoveredListNode))
			{
				listNode = this.hoveredListNode;
				htmlNode = this.getIconHtmlNodeFromListNode(listNode);
				event = "click";
			}
			else
			{
				listNode = this.getParentListNode(this.hoveredListNode)
				htmlNode = this.getHtmlNodeFromListNode(listNode);
				event = "mouseover";
			}
		}
		if(htmlNode != null)
		{
			if (e.preventDefault)
	            e.preventDefault();
			fireEvent(htmlNode, event);
		}
	}
	GloList.prototype.getHtmlNodeFromListNode = function(listNode)
	{
		if(listNode == null)
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
	GloList.prototype.getListNodeFromHtmlNode = function(htmlNode)
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
	GloList.prototype.getLoadingHtmlNode = function()
	{
		var node = document.getElementById(GLOLOADINGDIVID);
		return node;
	}	
	GloList.prototype.getHtmlChildNode = function(node, tag, type)
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
		return null;
	}
	GloList.prototype.getFirstListNode = function()
	{
		var listNodes = this.rootHtmlNode.getElementsByTagName("div");
		if(listNodes.length > 0)
			return listNodes[0];
	}
	GloList.prototype.getNextTreeItem = function(listNode)
	{
		var nextNode = null;
		if(listNode == null || listNode == 'undefined')
			return null;
		if(this.isBookOpenListNode(listNode))
			nextNode = this.getFirstChildNode(listNode);
		if(nextNode == null)
			nextNode = listNode.nextSibling; 
		if(nextNode == null)
		{
			var parentBookNode = listNode;

			while(nextNode == null)
			{
				parentBookNode = this.getParentListNode(parentBookNode);
				if(parentBookNode == null)
					break;
				nextNode = this.getNextSiblingNode(parentBookNode);
			}
		}
		return nextNode;	
	}
	GloList.prototype.getPreviousTreeItem = function(listNode)
	{
		var prevNode = this.getPreviousSiblingNode(listNode);
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
			prevNode = this.getParentListNode(listNode);
		return prevNode;	
	}
	GloList.prototype.getFirstChildNode = function(listNode)
	{
		if(listNode == null || listNode == 'undefined')
			return null;
		var bookChildsNode = null;
		if(this.isBookOpenListNode(listNode))
			bookChildsNode = listNode.getElementsByTagName("div")[1];
		if(bookChildsNode != null)
			return bookChildsNode.firstChild;
		else
			return null;
	}
	GloList.prototype.getLastChildNode = function(listNode)
	{
		var lastChild = null;
		if(listNode == null || listNode == 'undefined')
			return null;
		lastChild = this.getFirstChildNode(listNode);
		nextSibling = lastChild;
		while(nextSibling != null)
		{
			lastChild = nextSibling;
			nextSibling = this.getNextSiblingNode(nextSibling);
		}
		return lastChild;
	}
	GloList.prototype.getParentListNode = function(listNode)
	{
		if(listNode == null || listNode == 'undefined')
			return null;
		var bookChildsNode = listNode.parentNode;
		if(bookChildsNode != null && this.isNodeItemTypeThis(bookChildsNode,ITEMTYPEBOOKCHILDS) == false)
			return null;
		else
			return bookChildsNode.parentNode;
	}
	GloList.prototype.getNextSiblingNode = function(listNode)
	{
		if(listNode == null || listNode == 'undefined')
			return null;
		var nextSibling = listNode.nextSibling;
		if(nextSibling == null || nextSibling.tagName == null || nextSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return nextSibling;
	}
	GloList.prototype.getPreviousSiblingNode = function(listNode)
	{
		if(listNode == null || listNode == 'undefined')
			return null;
		var prevSibling = listNode.previousSibling;		
		if(prevSibling == null || prevSibling.tagName == null || prevSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return prevSibling;
	}
	GloList.prototype.isBookOpenListNode = function(listNode)
	{
		if(listNode == null)
			return false;
		var src = listNode.getAttribute(DATASRC);
		if(src == null || src == '')
			return false;
		else
		{
			var htmlNode = this.getHtmlNodeFromListNode(listNode);
			if(this.isNodeItemTypeThis(htmlNode,ITEMTYPEBOOKOPEN))
				return true;
			else
				return false;
		}
	}
	GloList.prototype.isBookClosedListNode = function(listNode)
	{
		if(listNode == null)
			return false;
		var src = listNode.getAttribute(DATASRC);
		if(src == null || src == '')
			return false;
		else
		{
			var htmlNode = this.getHtmlNodeFromListNode(listNode);
			if(this.isNodeItemTypeThis(htmlNode,ITEMTYPEBOOKCLOSED))
				return true;
			else
				return false;
		}
	}	
	GloList.prototype.isUrlNode = function(listNode)
	{
		var htmlNode = this.getHtmlNodeFromListNode(listNode);
		if(this.isNodeItemTypeThis(htmlNode, ITEMTYPEURL))
			return true;
		else
			return false;	
	}
	GloList.prototype.setLoadingDisplayInfo = function(iconClass, iconHtml, textClass, textString)
	{
		this.loadingIconClass = iconClass;
		this.loadingIconHtml = iconHtml;
		this.loadingTextClass = textClass;
		this.loadingText = textString;
	}
	GloList.prototype.insertLoadingMsg = function(htmlNode)
	{
		var divLoading = document.createElement('div');
		divLoading.className = TREEITEMCLASS;
		this.setNodeItemType(divLoading, ITEMTYPELOADING);
		divLoading.setAttribute("id", GLOLOADINGDIVID);
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
	GloList.prototype.removeLoadingMsg = function(htmlNode)
	{
		var divLoading = this.getLoadingHtmlNode();
		htmlNode.removeChild(divLoading);	
	}
	GloList.prototype.hoverNode = function(htmlNode, hoverClass)
	{
		if(this.hoveredListNode != null)
		{
			var htmNode = this.getHtmlNodeFromListNode(this.hoveredListNode);
			fireEvent(htmNode, 'mouseout');
		}
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		if(listNode != this.selectedListNode)
			htmlNode.className = hoverClass + " " + UNSELECTABLECLASS;
		this.hoveredListNode = listNode;
	}
	GloList.prototype.focusHoveredNode = function()
	{
		if(this.hoveredListNode == null)
			this.hoveredListNode = this.getFirstListNode();
		if(this.hoveredListNode != null)
		{
			var htmlNode = this.getHtmlNodeFromListNode(this.hoveredListNode);
			fireEvent(htmlNode, 'mouseover');
		}
	}
	GloList.prototype.blurHoveredNode = function()
	{
		if(this.hoveredListNode != null)
		{
			var htmlNode = this.getHtmlNodeFromListNode(this.hoveredListNode);
			fireEvent(htmlNode, 'mouseout');
		}
	}
	GloList.prototype.hoverOutNode = function(htmlNode, normalClass)
	{
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		if(listNode != this.selectedListNode)
		{
			htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
		}
	}
	GloList.prototype.toggleNode = function(htmlNode)
	{
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		var defNode = this.getHtmlChildNode(listNode, "div", ITEMTYPEDEF)
		if(defNode.style.display == "none")
			defNode.style.display = "block";
		else
			defNode.style.display = "none";
	}
	GloList.prototype.clickNode = function(htmlNode, clickClass, normalClass, url)
	{
		var listNode = this.getListNodeFromHtmlNode(htmlNode);
		if(listNode == this.selectedListNode)
			return;
		htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
		document.location = url;
	}
	GloList.prototype.isNodeItemTypeThis = function(node, type)
	{
		if(this.getNodeItemType(node) == type)
			return true;
		else
			return false;
	
	}
	GloList.prototype.getNodeItemType = function(node)
	{
		if(node != null && node != 'undefined')
			return node.getAttribute(DATAITEMTYPE);
		else
			return null;
	}
	GloList.prototype.setNodeItemType = function(node, type)
	{
		if(node != null && node != 'undefined')
			node.setAttribute(DATAITEMTYPE, type);	
	}
	GloList.prototype.getNodeTerm = function(node)
	{
		if(node != null && node != 'undefined')
			return node.getAttribute(DATATERM);
		else
			return null;
	}
	GloList.prototype.setNodeTerm = function(node, term)
	{
		if(node != null && node != 'undefined')
			node.setAttribute(DATATERM, term);	
	}
}

function onGloNodeHover(htmlNode, hoverClass)
{
	gGloList.hoverNode(htmlNode, hoverClass);
}
function onGloNodeHoverOut(htmlNode, normalClass)
{
	gGloList.hoverOutNode(htmlNode, normalClass);
}
function onGloNodeClick(htmlNode)
{
	gGloList.toggleNode(htmlNode);
}
function onGloListKeyPress(e)
{
	gGloList.pressKey(e);
}
function onGloListFocus()
{
	gGloList.focusHoveredNode();
}
function onGloListBlur()
{
	gGloList.blurHoveredNode();
}
function callbackGloRootFileLoaded(xmlDoc, arg) //Cannot use binding as IE9 does not support it
{
	gGloList.loadRootFiles(xmlDoc, arg);
}
function callbackGloChunkLoaded(xmlDoc, arg)
{
	gGloList.readChunk(xmlDoc, arg);
}
function onGloNodeExpand(e, plusHtmlNode)
{
	var evt = e || window.event;
	if(evt.preventDefault)
	{
		evt.preventDefault();  
	}
	else
	{
		evt.returnValue = false;  
		evt.cancelBubble=true;  
	}
	var listNode = gGloList.getListNodeFromHtmlNode(plusHtmlNode.parentNode);
	gGloList.expandListNode(listNode);
}
function onGloNodeCollapse(e, minusHtmlNode)
{
	var evt = e || window.event;
	if(evt.preventDefault)
	{
		evt.preventDefault();  
	}
	else
	{
		evt.returnValue = false;  
		evt.cancelBubble=true;  
	}
	var listNode = gGloList.getListNodeFromHtmlNode(minusHtmlNode.parentNode);
	gGloList.collapseListNode(listNode);
}
function filterGlo(e)
{
	if(e != null && e.type == 'submit')
		preventEvent(e);
	return gGloList.filterKeywords();
}
function callbackGloCSHModeRead(cshmode)
{
	if(cshmode == CSHMODE)
	{
		var filterBox = document.getElementById(gGloList.filterBoxId);
		if(filterBox != null)
			patchInputForSubmit(filterBox, function(){filterGlo(event);});
	}
}