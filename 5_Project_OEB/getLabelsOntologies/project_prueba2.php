<?php
require 'vendor/autoload.php';

$json = getListOntologyForForm();

echo $json;

//$error = json_last_error();

//var_dump($json, $error === JSON_ERROR_UTF8);

function getListOntologyForForm() {
	//variables
	$resource;
	$graph;
	$classArray = array();
    $label;
    $nameUrlOntology = "oeb_formats";

	if (isset($nameUrlOntology)) {

		//ontology-general
		$graph = EasyRdf_Graph::newAndLoad("https://raw.githubusercontent.com/inab/OEB-ontologies/master/oebDatasets-complete.owl","rdfxml");
		
		if ($nameUrlOntology == "oeb_datasets") {
			$resource = $graph->resource("https://w3id.org/oebDatasets/dataset");
		} elseif ($nameUrlOntology == "oeb_formats") {
			$resource = $graph->resource("https://w3id.org/oebDataFormats/FormatDatasets");
		} 
		
 		//get all the classes that are subclass of the uri 'https://w3id.org/oebDataFormats/FormatDatasets'
		$classes = $graph->resourcesMatching("rdfs:subClassOf",$resource);

		//get the first classes (without showing the childrens)
		foreach ($classes as $class) {
			//get the label of first classes (without showing the childrens)
            $label = $class->getLiteral('rdfs:label');
            $label = (string)$label; 

			//get the uri of classes that extends from the previous class find
			$resourceClassesInherited = $graph->resource($class);

			//get all the classes that are subclass of the uri found in the previous step (all the uris of classes that extends from the first classes found)
			$classesInherited = $graph->resourcesMatching("rdfs:subClassOf",$resourceClassesInherited);
			
			//if there are not any format inherited in the first classes do not do it
			if ($classesInherited != null) {
				//get the classes inherited (the childrens)
				foreach($classesInherited as $classInherited) {
					//get the label of the classes inherited (the childrens)
                    $labelClassInherited = $classInherited->getLiteral('rdfs:label');
                    $labelClassInherited = (string)$labelClassInherited;
					array_push($classArray, $labelClassInherited);
				}
			} else {
				array_push($classArray, $label);
			}
        }

        $array = array("JSON", "PDB");

        $array_def = array(
            "label" => $array
        );

        $process_json = json_encode($classArray);
        
		return $process_json;
	} else {
		return "nohola";
		//is not exist = not authorized
	}
}
