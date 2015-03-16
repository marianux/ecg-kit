
function IdxRootFileXmlObject(xmlDoc, i)
{
	this.xmlDoc = xmlDoc;
	this.nodeIndex = i;
}
function IdxChunkXmlObject(xmlDoc, i, len)
{
	this.xmlDoc = xmlDoc;
	this.nodeIndex = i;
	this.length = len;
}
function KeywordInfoObject(node, path, parentObj)
{
	this.node = node;
	this.path = path;
	this.nodeIndex = 0;
	this.parentObj = parentObj;
}
function IdxTree(idxRootPathsArr, dataFolder, rootFile)
{
	this.rootPathsArr = idxRootPathsArr;
	this.rootFilesXmlArr = new Array;
	this.curChunksXmlArr = new Array;
	this.nextChunkIndex = 0;
	this.curKeyName = "";
	this.curChunksToBeMerged = new Array;
	this.curCategory = "";

	this.errorMsg = "";
	this.kWClass = "";
	this.kWStyle = "";
	this.kWClassHover = "";
	this.kWClassClick = "";
	this.linkClass = "";
	this.linkStyle = "";
	this.linkClassHover = "";
	this.linkClassClick = "";
	this.categoryClass = "";
	this.categoryStyle = "";
	this.kWHtml = "";
	this.linkHtml = "";
	this.categoryHtml = "";
	
	this.bookChildsClass = "";
	this.filterBoxId = "";
	
	this.dataFolder = dataFolder;
	this.saveNodesState = true;
	this.rootHtmlNode = null;
	
	this.rootFile = rootFile;
	this.selectedTreeNode = null;
	this.hoveredTreeNode = null;
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
	
	IdxTree.prototype.init = function()
	{
		this.rootHtmlNode = document.getElementById(this.rootId);
		if(this.rootHtmlNode.attachEvent)
		{
			this.rootHtmlNode.attachEvent("onkeydown", function(){onIdxTreeKeyPress(event);});
			this.rootHtmlNode.attachEvent("onblur", function(){onIdxTreeBlur();});
			this.rootHtmlNode.attachEvent("onfocus", function(){onIdxTreeFocus();});
		}
		else
		{
			this.rootHtmlNode.setAttribute("onkeydown" , "onIdxTreeKeyPress(event)");	
			this.rootHtmlNode.setAttribute("onblur" , "onIdxTreeBlur()");	
			this.rootHtmlNode.setAttribute("onfocus" , "onIdxTreeFocus()");	
		}
		if(gbIE5)
			readSetting(RHCSHMODE, callbackIdxCSHModeRead);
	}
	IdxTree.prototype.load = function()
	{
		this.insertLoadingMsg(this.rootHtmlNode);
		var i=0;
		xmlJsReader.loadFile(this.rootPathsArr[0] + "/" + this.dataFolder + "/" + this.rootFile, callbackRootFileLoaded, i);
	}
	IdxTree.prototype.loadRootFiles = function(xmlDoc, i)
	{
		var rootFileXmlObj = new IdxRootFileXmlObject(xmlDoc, 0);
		this.rootFilesXmlArr[i] = rootFileXmlObj;
		chunkXmlObj = new IdxChunkXmlObject(null, 0, 0);
		this.curChunksXmlArr[i] = chunkXmlObj;
		++i;
		var len = this.rootPathsArr.length;
		if(i < len)
			xmlJsReader.loadFile(this.rootPathsArr[i] + "/" + this.dataFolder + "/" + this.rootFile, callbackRootFileLoaded, i);
		else
			this.mergeKeywords();
	}
	IdxTree.prototype.mergeKeywords = function()
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
			var indexXmlNode = rootFileXmlObj.xmlDoc.getElementsByTagName(INDEXNODE)[0];
			var childNodes = indexXmlNode.getElementsByTagName(CHUNKINFONODE);
			var len = childNodes.length;
			if(rootFileXmlObj.nodeIndex < len)
			{
				var chunkInfoNode = childNodes[rootFileXmlObj.nodeIndex];
				rootFileXmlObj.nodeIndex++;
				var url = chunkInfoNode.getAttribute(URL);
				var chunkPath = this.rootPathsArr[this.nextChunkIndex] + "/" + this.dataFolder + "/" + url;
				xmlJsReader.loadFile(chunkPath, callbackChunkLoaded);
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
			var childNodes = getChildElementsByTagName(dataNode,KEYNODE);
			var keyNode = childNodes[chunkXmlObj.nodeIndex];
			var keyName = keyNode.getAttribute(NAME);
			if(this.curKeyName == "" || compare(keyName,this.curKeyName) == -1)
			{
				this.curKeyName = keyName;
				this.curChunksToBeMerged.splice(0, this.curChunksToBeMerged.length);
				var kwInfoObj = new KeywordInfoObject(keyNode, this.rootPathsArr[this.nextChunkIndex], chunkXmlObj);
				this.curChunksToBeMerged[0] = kwInfoObj;
			}
			else if(keyName == this.curKeyName)
			{
				var kwInfoObj = new KeywordInfoObject(keyNode, this.rootPathsArr[this.nextChunkIndex], chunkXmlObj);
				this.curChunksToBeMerged[this.curChunksToBeMerged.length] = kwInfoObj;
			}
			this.incrementNextChunkIndexAndDoAction();
			return;
		}
	}
	IdxTree.prototype.incrementNextChunkIndexAndDoAction = function()
	{
		this.nextChunkIndex++;
		if(this.nextChunkIndex >= this.rootPathsArr.length)
		{
			this.nextChunkIndex = 0;
			if(this.curChunksToBeMerged.length != 0)
			{
				this.insertKeyword(this.rootHtmlNode, this.curChunksToBeMerged, ITEMTYPEKW);
				this.curChunksToBeMerged.splice(0, this.curChunksToBeMerged.length);
				this.curKeyName = "";
			}
			else
			{
				this.removeLoadingMsg(this.rootHtmlNode);
				this.curKeyName = "";
				this.filterKeywords(true);
				return;
			}
		}
		this.mergeKeywords();
	}
	IdxTree.prototype.readChunk = function(xmlDoc, arg)
	{
		var chunkXmlObj = this.curChunksXmlArr[this.nextChunkIndex];
		var dataNode = xmlDoc.getElementsByTagName(DATANODE)[0];
		chunkXmlObj.xmlDoc = xmlDoc;
		chunkXmlObj.nodeIndex = 0;
		chunkXmlObj.length = getChildElementsByTagName(dataNode,KEYNODE).length;
		this.mergeKeywords();
	}
	IdxTree.prototype.insertKeyword = function(parentHtmlNode, chunksArr, itemType)
	{
		var topicName = "";
		var topicNode = null;
		var path = null;
		var url = null;
		var keyNode = null;
		var kwName = "";
		
		
		var classNormal = this.kWClass;
		var classHover = this.kWClassHover;
		var classClick = this.kWClassClick;
		var inlinestyle = this.kWStyle;
		var classChilds = this.bookChildsClass;
		var linkClassNormal = this.linkClass;
		var linkClassHover = this.linkClassHover;
		var linkClassClick = this.linkClassClick;
		var linkInlineStyle = this.linkStyle;
		var html = this.kWHtml;
		var linkHtml = this.linkHtml;
		
		var len = chunksArr.length;
		
		for(var i=0; i<len; i++)
		{
			var kwInfoObj = chunksArr[i];
			if(kwInfoObj.parentObj != null)
				kwInfoObj.parentObj.nodeIndex++;
		}
		
		var noTopics = 0;
		var path = "";
		var topicNodes = null;
		var topicNode = null;
		for(var i=0; i<len; i++)
		{
			var kwInfoObj = chunksArr[i];
			keyNode = kwInfoObj.node;
			kwName = keyNode.getAttribute(NAME);
			path = kwInfoObj.path;
			var tempTopicsNodes = getChildElementsByTagName(keyNode,TOPICNODE);
			if(tempTopicsNodes.length > 0)
			{
				topicNodes = tempTopicsNodes;
				noTopics += topicNodes.length;
			}				
		}
		var ch = kwName.substring(0, 1).toLocaleUpperCase();
		if(itemType == ITEMTYPEKW && compare(ch,this.curCategory)!=0)
		{
			this.curCategory = ch;
			this.insertCategory(ch);
		}
		
		var treeNode = document.createElement("div");
		treeNode.setAttribute('class', TREEITEMCLASS);
		if(parentHtmlNode == this.rootHtmlNode)
		{
			var divLoading = this.getLoadingHtmlNode();
			parentHtmlNode.insertBefore(treeNode, divLoading);	
		}
		else
			parentHtmlNode.appendChild(treeNode);
		
		if(noTopics == 1)
		{
			var topicUrl = topicNodes[0].getAttribute(URL);
			if(topicUrl != null && !_isHTTPUrl(topicUrl))
				topicUrl = path + "/" + topicUrl;
			this.insertChildHtmlNode(treeNode, kwName, itemType, html, classNormal, classHover, classClick, inlinestyle, topicUrl);
		}	
		else
		{
			this.insertChildHtmlNode(treeNode, kwName, itemType, html, classNormal, classHover, classClick, inlinestyle, null);
			var curTopicName = "";
			var curPath = null;
			var curKwInfoObj = null;
	
			var childsHtmlNode = document.createElement("div");
			childsHtmlNode.className = classChilds;
			childsHtmlNode.style.display = "none";
			treeNode.appendChild(childsHtmlNode);
			this.setNodeItemType(childsHtmlNode, ITEMTYPEBOOKCHILDS);
			
			while(true)
			{
				curTopicName = "";
				curPath = null;
				curUrl = "";
				curKwInfoObj = null;
				for(var i=0; i<len; i++)
				{
					var kwInfoObj = chunksArr[i];
					if(kwInfoObj.nodeIndex == -1)
						continue;
					path = kwInfoObj.path;
					topicName = "";
					keyNode = kwInfoObj.node;
					topicNodes = getChildElementsByTagName(keyNode,TOPICNODE);
					if(kwInfoObj.nodeIndex < topicNodes.length)
					{
						topicNode = topicNodes[kwInfoObj.nodeIndex];
						topicName = topicNode.getAttribute(NAME);
						url = topicNode.getAttribute(URL);
					}
					else
					{
						kwInfoObj.nodeIndex = -1;
						continue;
					}
					if(curTopicName == "" || compare(topicName, curTopicName) == -1)
					{
						curTopicName = topicName;
						curPath = path;
						curUrl = url;
						curKwInfoObj = kwInfoObj;
					}
				}
				if(curTopicName == "" && curPath == null && curUrl == "" && curKwInfoObj == null)
					break;
				curKwInfoObj.nodeIndex++;
				var topicUrl = null;
				if(url != null)
				{
					if(curUrl != null && !_isHTTPUrl(curUrl))
						topicUrl = curPath + "/" + curUrl;
					else
						topicUrl = curUrl;
				}
				this.insertChildHtmlNode(childsHtmlNode, curTopicName, ITEMTYPELINK, linkHtml, linkClassNormal, linkClassHover, linkClassClick, linkInlineStyle, topicUrl); 
			}
		}	
		var subKwParentHtmlNode = null;
		var subKwNode = null;
		var curSubKwsToBeMergedArr = new Array;
		var curSubKwName = "";
		while(true)
		{
			curSubKwsToBeMergedArr.splice(0, curSubKwsToBeMergedArr.length);
			curSubKwName = ""; 
			for(var i=0; i<len; i++)
			{
				var kwInfoObj = chunksArr[i];
				var node = kwInfoObj.node;
				if(node == null)
					continue;
				if(kwInfoObj.nodeIndex == -1)
					kwInfoObj.nodeIndex = 0;
				path = kwInfoObj.path;
				var subKwNodes = getChildElementsByTagName(node,KEYNODE);
				var subkwName = "";
				if(kwInfoObj.nodeIndex < subKwNodes.length)
				{
					subKwNode = subKwNodes[kwInfoObj.nodeIndex];
					subkwName = subKwNode.getAttribute(NAME);
				}
				else
				{
					kwInfoObj.node = null;
					continue;
				}
				if(curSubKwName == "" || compare(subkwName, curSubKwName) == -1)
				{
					curSubKwName = subkwName;
					curSubKwsToBeMergedArr.splice(0, curSubKwsToBeMergedArr.length);
					var subKwInfoObj = new KeywordInfoObject(subKwNode, path, kwInfoObj);
					curSubKwsToBeMergedArr[0] = subKwInfoObj;
				}
				else if(subkwName == curSubKwName)
				{
					var subKwInfoObj = new KeywordInfoObject(subKwNode, path, kwInfoObj);
					curSubKwsToBeMergedArr[curSubKwsToBeMergedArr.length] = subKwInfoObj;
				}
			}
			if(curSubKwsToBeMergedArr.length > 0)
			{
				if(subKwParentHtmlNode == null)
				{
					subKwParentHtmlNode = document.createElement("div");
					subKwParentHtmlNode.className = classChilds;
					treeNode.appendChild(subKwParentHtmlNode);
				}
				this.insertKeyword(subKwParentHtmlNode, curSubKwsToBeMergedArr, ITEMTYPESUBKW);
			}
			else
				break;
		}
	}
	IdxTree.prototype.insertChildHtmlNode = function(parentHtmlNode, name, itemType, html, classNormal, classHover, classClick, style, url)
	{
		var bAddAnchor = false;
		if(url != null && url != "")
			bAddAnchor = true;
		html = html.replace(LINK_NAME_MACRO, name);

		var htmlNode = document.createElement("div");
		htmlNode.className = classNormal + " " + UNSELECTABLECLASS;
		if(style != "")
			htmlNode.style.cssText = style;
		htmlNode.setAttribute("title", name);
		htmlNode.innerHTML = html;
		if(bAddAnchor)
		{
			var anchorNode = document.createElement("a");
			anchorNode.className = NOLINKANCHORCLASS;
			anchorNode.appendChild(htmlNode);
			parentHtmlNode.appendChild(anchorNode);
		}
		else
			parentHtmlNode.appendChild(htmlNode);
		var urlWithId = this.addEventsToNode(htmlNode, classNormal, classHover, classClick, url);
		if(bAddAnchor)
			anchorNode.setAttribute("href", url);	
		this.setNodeItemType(htmlNode, itemType);
		if(itemType == ITEMTYPEKW)
			this.setNodeTerm(htmlNode, name);
	}
	IdxTree.prototype.insertCategory = function(ch)
	{
		var treeNode = document.createElement("div");
		treeNode.setAttribute('class', TREEITEMCLASS);
		
		var categoryElem = document.createElement("div");
		categoryElem.className = this.categoryClass;
		categoryElem.style.cssText = this.categoryStyle;
		categoryElem.innerHTML = this.categoryHtml.replace(LINK_NAME_MACRO, ch);
		this.setNodeItemType(categoryElem, ITEMTYPECATEGORY);
		treeNode.appendChild(categoryElem);
		var divLoading = this.getLoadingHtmlNode();
		this.rootHtmlNode.insertBefore(treeNode, divLoading);
	}
	IdxTree.prototype.addEventsToNode = function(htmlNode, classNormal, classHover, classClick, url)
	{
		if(htmlNode.attachEvent)
		{
			if(isTouchDevice())
			{
				htmlNode.attachEvent('ontouchstart', function(){onIdxNodeHover(htmlNode, classHover);});
				htmlNode.attachEvent('ontouchend', function(){onIdxNodeHoverOut(htmlNode, classHover);});
				htmlNode.attachEvent('ontouchmove', function(){onIdxNodeHoverOut(htmlNode, classHover);});
			}
			else
			{
				htmlNode.attachEvent('onmouseout', function(){onIdxNodeHoverOut(htmlNode, classNormal);});
				htmlNode.attachEvent('onmouseover', function(){onIdxNodeHover(htmlNode, classHover);});
			}
			if(url == null)
				htmlNode.attachEvent('onclick', function(){onIdxNodeClick(htmlNode);});
		}
		else
		{
			if(isTouchDevice())
			{
				htmlNode.setAttribute("ontouchstart", "onIdxNodeHover(this,'" + classHover + "')");
				htmlNode.setAttribute("ontouchend", "onIdxNodeHoverOut(this,'" + classNormal + "')");
			}
			else
			{
				htmlNode.setAttribute("onmouseout", "onIdxNodeHoverOut(this,'" + classNormal + "')");
				htmlNode.setAttribute("onmouseover", "onIdxNodeHover(this,'" + classHover + "')");
			}
			if(url == null)
				htmlNode.setAttribute("onclick", "onIdxNodeClick(this)");
		}
	}
	IdxTree.prototype.filterKeywords = function(bIsOnLoad)
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
			treeNode = listNodes[i];
			if(treeNode.nodeType != JS_TAGTOKEN)
				continue;
			var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
			var itemType = this.getNodeItemType(htmlNode);
			if(itemType == ITEMTYPEKW)
			{
				var strItemValue = this.getNodeTerm(htmlNode).toLocaleLowerCase();
				if(strItemValue.indexOf(strFilter)>=0)
				{
					treeNode.style.display = "block";
					foundDisplayNode = true;
					itemDisplayed = true;
				}
				else
					treeNode.style.display = "none";
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
				categoryNode = treeNode;
				foundDisplayNode = false;
			}
			else if(itemType == ITEMTYPEBOOKCHILDS)
			{
				if(itemDisplayed == true)
					treeNode.style.display = "block";
				else
					treeNode.style.display = "none";
				itemDisplayed = false;
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
	IdxTree.prototype.pressKey = function(e)
	{

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
		{
			if (e.preventDefault)
	            e.preventDefault();
			fireEvent(htmlNode, event);
		}
	}
	IdxTree.prototype.getHtmlNodeFromTreeNode = function(treeNode)
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
	IdxTree.prototype.getTreeNodeFromHtmlNode = function(htmlNode)
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
	IdxTree.prototype.getLoadingHtmlNode = function()
	{
		var node = document.getElementById(IDXLOADINGDIVID);
		return node;
	}
	IdxTree.prototype.getHtmlChildNode = function(node, tag, type)
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
	IdxTree.prototype.getFirstTreeNode = function()
	{
		var treeNodes = this.rootHtmlNode.getElementsByTagName("div");
		if(treeNodes.length > 0)
			return treeNodes[0];
	}
	IdxTree.prototype.getNextTreeItem = function(treeNode)
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
	IdxTree.prototype.getPreviousTreeItem = function(treeNode)
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
	IdxTree.prototype.getFirstChildNode = function(treeNode)
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
	IdxTree.prototype.getLastChildNode = function(treeNode)
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
	IdxTree.prototype.getParentTreeNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var bookChildsNode = treeNode.parentNode;
		if(bookChildsNode != null && this.isNodeItemTypeThis(bookChildsNode,ITEMTYPEBOOKCHILDS) == false)
			return null;
		else
			return bookChildsNode.parentNode;
	}
	IdxTree.prototype.getNextSiblingNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var nextSibling = treeNode.nextSibling;
		if(nextSibling == null || nextSibling.tagName == null || nextSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return nextSibling;
	}
	IdxTree.prototype.getPreviousSiblingNode = function(treeNode)
	{
		if(treeNode == null || treeNode == 'undefined')
			return null;
		var prevSibling = treeNode.previousSibling;		
		if(prevSibling == null || prevSibling.tagName == null || prevSibling.tagName.toLowerCase() != 'div')
			return null;
		else
			return prevSibling;
	}
	IdxTree.prototype.isBookOpenTreeNode = function(treeNode)
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
	IdxTree.prototype.isBookClosedTreeNode = function(treeNode)
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
	IdxTree.prototype.isUrlNode = function(treeNode)
	{
		var htmlNode = this.getHtmlNodeFromTreeNode(treeNode);
		if(this.isNodeItemTypeThis(htmlNode, ITEMTYPEURL))
			return true;
		else
			return false;	
	}
	IdxTree.prototype.setLoadingDisplayInfo = function(iconClass, iconHtml, textClass, textString)
	{
		this.loadingIconClass = iconClass;
		this.loadingIconHtml = iconHtml;
		this.loadingTextClass = textClass;
		this.loadingText = textString;
	}
	IdxTree.prototype.insertLoadingMsg = function(htmlNode)
	{
		var divLoading = document.createElement('div');
		divLoading.className = TREEITEMCLASS;
		divLoading.setAttribute("id", IDXLOADINGDIVID);
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
	IdxTree.prototype.removeLoadingMsg = function(htmlNode)
	{
		var divLoading = this.getLoadingHtmlNode();
		htmlNode.removeChild(divLoading);	
	}
	IdxTree.prototype.hoverNode = function(htmlNode, hoverClass)
	{
		if(this.hoveredTreeNode != null)
		{
			var htmNode = this.getHtmlNodeFromTreeNode(this.hoveredTreeNode);
			fireEvent(htmNode, 'mouseout');
		}
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		if(treeNode != this.selectedTreeNode)
			htmlNode.className = hoverClass + " " + UNSELECTABLECLASS;
		this.hoveredTreeNode = treeNode;
	}
	IdxTree.prototype.focusHoveredNode = function()
	{
		if(this.hoveredTreeNode == null)
			this.hoveredTreeNode = this.getFirstTreeNode();
		if(this.hoveredTreeNode != null)
		{
			var htmlNode = this.getHtmlNodeFromTreeNode(this.hoveredTreeNode);
			fireEvent(htmlNode, 'mouseover');
		}
	}
	IdxTree.prototype.blurHoveredNode = function()
	{
		if(this.hoveredTreeNode != null)
		{
			var htmlNode = this.getHtmlNodeFromTreeNode(this.hoveredTreeNode);
			fireEvent(htmlNode, 'mouseout');
		}
	}
	IdxTree.prototype.hoverOutNode = function(htmlNode, normalClass)
	{
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		if(treeNode != this.selectedTreeNode)
		{
			htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
		}
	}
	IdxTree.prototype.toggleNode = function(htmlNode)
	{
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		var bookChildsNode = this.getHtmlChildNode(treeNode, "div", ITEMTYPEBOOKCHILDS);
		if(bookChildsNode.style.display == "none")
			bookChildsNode.style.display = "block";
		else
			bookChildsNode.style.display = "none";
	}
	IdxTree.prototype.clickNode = function(htmlNode, clickClass, normalClass, url)
	{
		var treeNode = this.getTreeNodeFromHtmlNode(htmlNode);
		if(treeNode == this.selectedTreeNode)
			return;
		htmlNode.className = normalClass + " " + UNSELECTABLECLASS;
		document.location = url;
	}
	IdxTree.prototype.isNodeItemTypeThis = function(node, type)
	{
		if(this.getNodeItemType(node) == type)
			return true;
		else
			return false;
	
	}
	IdxTree.prototype.getNodeItemType = function(node)
	{
		if(node != null && node != 'undefined')
			return node.getAttribute(DATAITEMTYPE);
		else
			return null;
	}
	IdxTree.prototype.setNodeItemType = function(node, type)
	{
		if(node != null && node != 'undefined')
			node.setAttribute(DATAITEMTYPE, type);	
	}
	IdxTree.prototype.getNodeTerm = function(node)
	{
		if(node != null && node != 'undefined')
			return node.getAttribute(DATATERM);
		else
			return null;
	}
	IdxTree.prototype.setNodeTerm = function(node, term)
	{
		if(node != null && node != 'undefined')
			node.setAttribute(DATATERM, term);	
	}
}

function onIdxNodeHover(htmlNode, hoverClass)
{
	gIdxTree.hoverNode(htmlNode, hoverClass);
}
function onIdxNodeHoverOut(htmlNode, normalClass)
{
	gIdxTree.hoverOutNode(htmlNode, normalClass);
}
function onIdxNodeClick(htmlNode)
{
	gIdxTree.toggleNode(htmlNode);
}
function onIdxTreeKeyPress(e)
{
	gIdxTree.pressKey(e);
}
function onIdxTreeFocus()
{
	gIdxTree.focusHoveredNode();
}
function onIdxTreeBlur()
{
	gIdxTree.blurHoveredNode();
}
function callbackRootFileLoaded(xmlDoc, arg) //Cannot use binding as IE9 does not support it
{
	gIdxTree.loadRootFiles(xmlDoc, arg);
}
function callbackChunkLoaded(xmlDoc, arg)
{
	gIdxTree.readChunk(xmlDoc, arg);
}
function filterIdx(e)
{
	if(e != null && e.type == 'submit')
		preventEvent(e);
	return gIdxTree.filterKeywords();
}
function callbackIdxCSHModeRead(cshmode)
{
	if(cshmode == CSHMODE)
	{
		var filterBox = document.getElementById(gIdxTree.filterBoxId);
		if(filterBox != null)
			patchInputForSubmit(filterBox, function(){filterIdx(event);});
	}
}