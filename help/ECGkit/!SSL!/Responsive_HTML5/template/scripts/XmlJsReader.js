var gXMLBuffer = null;
var gFileNameToXMLMap = new Object();
var xmlJsReader = new XmlJsReader();

function XmlInfo(xmlPath, oFunCallback, args) {
    this.sXmlPath = xmlPath;
    this.callback = oFunCallback;
    this.cbargs = args;
}

function XmlJsReader() {
    this.queue = new MhQueue();
    this.bLoading = false;

    this.getJsNameFromXmlName = function (xmlPath) {
        var indx = xmlPath.lastIndexOf(".xml");
        if (indx != -1) {
            var jsPath = xmlPath.substring(0, indx);
            jsPath += "_xml.js";
            return jsPath;
        }
		return xmlPath;
    }
    /*use relative path for xmlPath*/
    this.loadFile = function (xmlPath, oFunCallback, args) {
        this.queue.enqueue(new XmlInfo(xmlPath, oFunCallback, args));
        this.loadFromQueue();
    }

    this.loadFromQueue = function () {
        if (this.queue.isEmpty() || this.bLoading) {
            return;
        }
        else {
            var xmlInfo = this.queue.peek();
            if (typeof (gFileNameToXMLMap[xmlInfo.sXmlPath]) == 'undefined') {
                var jsPath = this.getJsNameFromXmlName(xmlInfo.sXmlPath);
                this.loadScript(jsPath, this.onScriptLoaded);
            }
            else {
                this.onScriptLoaded();
            }
        }
    }

    this.onScriptLoaded = function () {
        var xmlInfo = xmlJsReader.queue.dequeue();
        if (typeof(gFileNameToXMLMap[xmlInfo.sXmlPath]) == 'undefined' && gXMLBuffer != null) {
            gFileNameToXMLMap[xmlInfo.sXmlPath] = gXMLBuffer;
        }
        var xmlDoc = null;
        if (typeof (gFileNameToXMLMap[xmlInfo.sXmlPath]) != 'undefined') {
            if (window.DOMParser) {
                var parser = new DOMParser();
                xmlDoc = parser.parseFromString(gFileNameToXMLMap[xmlInfo.sXmlPath], "text/xml");
            }
            else {
                xmlDoc = new ActiveXObject("Microsoft.XMLDOM");
                xmlDoc.async = false;
                var indx = gFileNameToXMLMap[xmlInfo.sXmlPath].indexOf("<?xml");
                if (indx != -1) {
                    indx = gFileNameToXMLMap[xmlInfo.sXmlPath].indexOf("?>", indx);
                    if (indx != -1) {
                        var strXML = gFileNameToXMLMap[xmlInfo.sXmlPath].substr(indx + 2);
                        xmlDoc.loadXML(strXML);
                    }
                }
                else {
                    xmlDoc.loadXML(gFileNameToXMLMap[xmlInfo.sXmlPath]);
                }
            }
        }
        gXMLBuffer = null;
        xmlJsReader.bLoading = false;

        if (xmlInfo.callback)
            xmlInfo.callback(xmlDoc, xmlInfo.cbargs);

        xmlJsReader.loadFromQueue();
    }

    this.loadScript = function (sScriptSrc, onScriptLoadedCB) {
        this.bLoading = true;
        var oHead = document.getElementsByTagName('head')[0];
        var oScript = document.createElement('script');
        oScript.type = 'text/javascript';
        oScript.charset = "utf-8";
		oScript.src = sScriptSrc;

        // IE 6 & 7
        if (oScript.readyState) {
            oScript.onreadystatechange = function () {
                if (oScript.readyState == 'loaded' ||
                    oScript.readyState == 'complete') {
                    onScriptLoadedCB();
                }
            }
        }
        else {
            oScript.onload = onScriptLoadedCB;
            oScript.onerror = onScriptLoadedCB;
        }

        oHead.appendChild(oScript);
    }
}
