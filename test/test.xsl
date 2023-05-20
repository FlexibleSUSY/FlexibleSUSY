<?xml version="1.0"?>
<xsl:stylesheet
    version="1.0" 
    xmlns="http://www.w3.org/1999/xhtml"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="xml" version="1.0" encoding="utf-8"
	      doctype-public="-//W3C//DTD XHTML 1.1//EN"
	      doctype-system="http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd"
	      omit-xml-declaration="no"
	      indent="yes"/>

  <xsl:template match="/">
    <html version="-//W3C//DTD XHTML 1.1//EN"
	  xml:lang="en"
	  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	  xsi:schemaLocation="http://www.w3.org/1999/xhtml
                              http://www.w3.org/MarkUp/SCHEMA/xhtml11.xsd"
	  >
      <head>
	<title>FlexibleSUSY test results</title>
	<meta name="author" content="Author"/>
	<meta name="description" content="FlexibleSUSY test results"/>
	<meta name="keywords" content="FlexibleSUSY"/>
	<meta http-equiv="content-type" content="application/xhtml+xml; charset=UTF-8"/>
        <style>
          table {
            text-align: left;
            margin-bottom: 1em;
            margin-top: 0.3em;
            border-collapse: collapse;
          }
          table {
             border-bottom: 2px solid #000;
          }
          table th {
             border-bottom: 1px solid #000;
             border-top: 2px solid #000;
          }
          td {
             padding: 0.2em 1em;
             margin: 5px 0 0;
          }
          th {
             padding: 0.4em 1em;
          }
          #ok {
             background-color: #0f0;
          }
          #failed {
             background-color: #f00;
          }
          a {
             color: #0022ff;
          }
          a:hover, a:focus, a:active {
             background-color: #0044ff;
             color: #fff;
             text-decoration: none;
          }
          tr:hover td {
             background-color: #ccc !important;
             color: #000;
          }
        </style>
      </head>
      <body>
	<h1>FlexibleSUSY test results</h1>
        <p><xsl:value-of select="/tests/@date"/></p>
        <table>
          <tr>
            <th>Test</th>
            <th>Status</th>
            <th>Date</th>
            <th>Commit</th>
          </tr>

          <xsl:for-each select="document(/tests/test/@filename)/test">
            <xsl:sort
               select="name"
               data-type="text"
               order="ascending"/>
            <tr>
              <td>
                <xsl:element name="a">
                  <xsl:attribute name="href">
                    <xsl:value-of select="logfile"/>
                  </xsl:attribute>
                  <xsl:value-of select="name"/>
                </xsl:element>
              </td>
	      <xsl:choose>
		<xsl:when test="status='0'">
                 <td id="ok">OK</td>
                </xsl:when>
		<xsl:otherwise>
                  <td id="failed">FAILED</td>
                </xsl:otherwise>
	      </xsl:choose>
              <td><xsl:value-of select="date"/></td>
              <td><xsl:value-of select="commit"/></td>
            </tr>
	  </xsl:for-each>

        </table>
      </body>
    </html>

  </xsl:template>

</xsl:stylesheet>
