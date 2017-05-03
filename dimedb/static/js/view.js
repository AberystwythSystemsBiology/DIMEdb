function basic_information(metabolite) {
    function fill_synonyms(synonyms_array) {
        if (synonyms_array != null) {
            var synonyms = "";
            for (i = 0; i < synonyms_array.length; i++) {
                synonyms += synonyms_array[i] + "; ";
            }
        }
        else {
            var synonyms = "Not available";
        }
        return synonyms
    }


    function fill_origins(origins_array) {
        if (origins_array.length != 0) {
            var origins = "";
            for (i = 0; i < origins_array.length; i++) {
                origins += origins_array[i] + "; ";
            }
        }

        else {
            var origins = "Not available";
        }

        return origins;
    }

    function fill_biofluids(biofluids_array) {
        if (biofluids_array != null) {
            var biofluids = "";
            for (i = 0; i < biofluids_array.length; i++) {
                biofluids += biofluids_array[i] + "; ";
            }
        }

        else {
            var biofluids = "Not available";
        }

        return biofluids;
    }

    function fill_tissue_locations(tissue_locations_array) {
        if (tissue_locations_array != null) {
            var tissues = "";
            for (i = 0; i < tissue_locations_array.length; i++) {
                tissues += tissue_locations_array[i] + "; ";
            }
        }
        else {
            var tissues = "Not available";
        }

        return tissues
    }

    $("#synonyms").html(fill_synonyms(metabolite["synonyms"]));
    $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
    $("#num_atoms").html(metabolite["num_atoms"]);
    $("#accurate_mass").html(metabolite["accurate_mass"].toFixed(6));
    $("#neutral_mass").html(metabolite["adducts"]["neutral"][0]["accurate_mass"].toFixed(6));
    $("#origins").html(fill_origins(metabolite["origins"]));
    $("#biofluids").html(fill_biofluids(metabolite["biofluid_location"]));
    $("#tissue").html(fill_tissue_locations(metabolite["tissue_locations"]))
}

function sidebar(metabolite) {
    function fill_pathways(pathways) {
        if (pathways != null) {
            for (pathway in pathways) {
                $("#kegg_pathways_list").append(
                    "<li class='list-group-item'>"+pathways[pathway]+"</li>"
                )
            }
        }

        else {
            $("#kegg_pathways_list").append(
                    "<li class='list-group-item'>None available.</li>"
                )

        }
    }

    function fill_sources(sources) {
        for (source in sources) {

            if (source == "kegg_id") {
                var source_url = "http://www.genome.jp/dbget-bin/www_bget?" + sources[source];
                var source_name = "KEGG";
                var source_colour = "list-group-item-success";
            }

            else if (source == "chebi_id") {
                var source_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:" + sources[source];
                var source_name = "CHEBI";
                var source_colour = "list-group-item-info";
            }

            else if (source == "hmdb_id") {
                var source_url = "http://www.hmdb.ca/metabolites/" + sources[source];
                var source_name = "HMDB";
                var source_colour = "list-group-item-warning";
            }

            else if (source == "pubchem_id") {
                var source_url = "https://pubchem.ncbi.nlm.nih.gov/compound/" + sources[source];
                var source_name = "PubChem";
                var source_colour = "list-group-item-danger";
            }

            var source_id = sources[source]

            $("#data_sources").append(
                "<li><b>"+source_name+"</b> (<a href='"+source_url+"' target='_blank'>"+source_id+"</a>)</li>"
            );
        }
    }
    fill_pathways(metabolite["pathways"])
    fill_sources(metabolite["sources"])
}

function isotopic_distributions(adducts) {
    function fill_adduct_panels(adducts) {
        var ionisations = ["negative", "positive"];

        for (i in ionisations) {
            var ionisation = ionisations[i];
            for (j in adducts[ionisation]) {
                var adduct  = adducts[ionisation][j];

                var button = "<div class='btn btn-primary btn-sm pull-right'>View Isotopes <i class='glyphicon glyphicon-signal'></i></div><div class='clearfix'></div>";

                $("#"+ionisation+"_adducts_list").append(
                    "<li class='list-group-item'><b>"+ adduct["type"] + ":</b> "+adduct["accurate_mass"].toFixed(4)+ button +"</li>"
                )
            }
        }

    }

    fill_adduct_panels(adducts);

}

function render_adduct_chart_and_table(adducts, ionisation, type) {

    function render_title(ionisation, type) {
        $("#adduct_label").html(ionisation + " " + type);
    }

    function render_distribution_table(i_d) {
        $("#distribution_table tbody tr").remove();
        for (i in i_d) {
            $("#distribution_table").append("<tr><td class='text-center'>"+i_d[i][0].toFixed(4)+"</td><td class='text-center'>"+i_d[i][1].toFixed(3)+"</td></tr>")
        }
    }


    function render_chart(distributions) {

        var x_plot = [];
        var y_plot = [];

        for (i in distributions) {
            var spectra = distributions[i];
            x_plot.push(spectra[0]);
            y_plot.push(spectra[1]);
        }

        var chart_data = [{
            x: x_plot,
            y: y_plot,
            type: 'bar',
            marker: {
                color: 'rgba(0, 0, 0, 1)'
            }
        }];

        var layout = {
            xaxis: {
                title: 'Mass-to-ion (m/z)',
                showgrid: false,
                range: [Math.min.apply(Math, x_plot) - 0.5, Math.max.apply(Math, x_plot) + 1]
            },
            yaxis: {
                title: 'Relative Intensity (%)'
            },
            margin: {
                l: 50,
                r: 50,
                b: 50,
                t: 50,
                pad: 1
            },
            bargap: 0.99
        };

        Plotly.newPlot("isotopic_distribution_chart", chart_data, layout, {displayModeBar: false});
    }

    render_title(ionisation, type);


    for (i in adducts[ionisation]) {
        if (adducts[ionisation][i]["type"] == type) {
            var adduct = adducts[ionisation][i];
            render_distribution_table(adduct["isotopic_distribution"])
            render_chart(adduct["isotopic_distribution"])
        }
    }
}

function render_metabolite_view(metabolite_id) {
    var metabolite = get_metabolite(metabolite_id);

    $("#name").html(metabolite["name"]);

    basic_information(metabolite);
    sidebar(metabolite);
    render_adduct_chart_and_table(metabolite["adducts"], "neutral", "[M]");
    isotopic_distributions(metabolite["adducts"]);

}