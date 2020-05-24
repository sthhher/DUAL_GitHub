<?php
require 'vendor/autoload.php';


//ontology-general
$graph = EasyRdf_Graph::newAndLoad("https://raw.githubusercontent.com/inab/OEB-ontologies/master/oebDatasets-complete.owl","rdfxml");

//file_types Ex. FASTA -> get the uri of formats
$resourceFormats = $graph->resource("https://w3id.org/oebDataFormats/FormatDatasets");
//data_types Ex. aggregation -> get the uri of datasets
$resourceDatasets = $graph->resource("https://w3id.org/oebDatasets/dataset");

//get all the formats that are subclass of the uri 'https://w3id.org/oebDataFormats/FormatDatasets'
$formats = $graph->resourcesMatching("rdfs:subClassOf",$resourceFormats);
//get all the datasets that are subclass of the uri 'https://w3id.org/oebDatasets/dataset'
$datasets = $graph->resourcesMatching("rdfs:subClassOf",$resourceDatasets);

?>

<h1>OWL FILE_TYPE</h1>
<?php
    //get the first formats (without showing the childrens)
    foreach ($formats as $format) {
        //get the label of first formats (without showing the childrens)
        $labelFormat = $format->getLiteral('rdfs:label');
        //print all the first formats (without showing childrens)
        echo "<br><br><b>FIRST SUBCLASS:</b> <br>";
        echo "<li>".$format." -- ".$labelFormat."</li>\n";

        //get the uri of formats that extends from the previous formats find
        $resourceFormatsInherited = $graph->resource($format);

        //get all the formats that are subclass of the uri found in the previous step (all the uris of formats that extends from the first formats found)
        $formatsInherited = $graph->resourcesMatching("rdfs:subClassOf",$resourceFormatsInherited);
        
        //if there are not any format inherited in the first formats do not do it
        if ($formatsInherited != null) {
            echo "<b>SECOND SUBCLASS:</b> <br>";
            //get the formats inherited (the childrens)
            foreach($formatsInherited as $formatInherited) {
                //get the label of the formats inherited (the childrens)
                $labelFormatInherited = $formatInherited->getLiteral('rdfs:label');
                //print all the formats inherited (the childrens)
                echo "<li>".$formatInherited." -- ".$labelFormatInherited."</li>\n";
            }
        }
    };
?>
<br><br>

<h1>OWL DATA_TYPE</h1>
<?php
    //get the first formats (without showing the childrens)
    foreach ($datasets as $dataset) {
        //get the label of first formats (without showing the childrens)
        $labelDataset = $dataset->getLiteral('rdfs:label');
        //print all the first formats (without showing childrens)
        echo "<br><br><b>FIRST SUBCLASS:</b> <br>";
        echo "<li>".$dataset." -- ".$labelDataset."</li>\n";

        //get the uri of formats that extends from the previous formats find
        $resourceDatasetsInherited = $graph->resource($dataset);

        //get all the formats that are subclass of the uri found in the previous step (all the uris of formats that extends from the first formats found)
        $datasetsInherited = $graph->resourcesMatching("rdfs:subClassOf",$resourceDatasetsInherited);
        
        //if there are not any format inherited in the first formats do not do it
        if ($datasetsInherited != null) {
            echo "<b>SECOND SUBCLASS:</b> <br>";
            //get the formats inherited (the childrens)
            foreach($datasetsInherited as $datasetInherited) {
                //get the label of the formats inherited (the childrens)
                $labelDatasetInherited = $datasetInherited->getLiteral('rdfs:label');
                //print all the formats inherited (the childrens)
                echo "<li>".$datasetInherited." -- ".$labelDatasetInherited."</li>\n";
            }
        }
    };
?>