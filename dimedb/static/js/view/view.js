function array_to_text(a) {
        if (a.length != 0) {
            var text = "";
            for (i = 0; i < a.length; i++) {
                text += a[i] + "; ";
            }
        }
        else {
            var text = "<i>Not available.</i>";
        }
        return text
}

function text_to_none(t) {
    if (t != null) {
        return t
    }
    else {
        return "<i>Not available.</i>"
    }
}

function fill_identification_infomation(identification_infomation, metabolite_id) {
    document.title = identification_infomation["Name"] +" : DIMEdb - Direct Infusion MEtabolite Database";
    $("#name").html(identification_infomation["Name"]);
    $("#synonyms").html(array_to_text(identification_infomation["Synonyms"]));
    $("#systematic_name").html(text_to_none(identification_infomation["Systematic Name"]));
    $("#iupac_name").html(text_to_none(identification_infomation["IUPAC Name"]));
    $("#inchi").html(text_to_none(identification_infomation["InChI"]));
    $("#inchi_key").html(text_to_none(metabolite_id));
    $("#smiles").html(text_to_none(identification_infomation["SMILES"]));
    $("#molecular_formula").html(identification_infomation["Molecular Formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
}

function fill_physicochemical_properties(physicochemical_properties) {
    $("#molecular_weight").html(physicochemical_properties["Molecular Weight"].toFixed(3));
    $("#hbd").html(physicochemical_properties["Hydrogen Bond Donors"]);
    $("#hba").html(physicochemical_properties["Hydrogen Bond Acceptors"]);
    $("#rotatable_bonds").html(physicochemical_properties["Rotatable Bonds"]);
    $("#formal_charge").html(physicochemical_properties["Formal Charge"]);
    $("#fraction_of_sp3").html(physicochemical_properties["Fraction of SP3 Carbon"].toFixed(4));
    $("#clogp").html(physicochemical_properties["clogP"].toFixed(4));
    $("#polar_surface_area").html(physicochemical_properties["Polar Surface Area"].toFixed(4));
    $("#mr_values").html(physicochemical_properties["MR Values"].toFixed(4));
    $("#rings").html(physicochemical_properties["Rings"]);
    $("#aromatic_rings").html(physicochemical_properties["Aromatic Rings"]);
    $("#heavy_atoms").html(physicochemical_properties["Heavy Atoms"]);
    $("#secondary_amines").html(physicochemical_properties["Secondary Amines"]);
    $("#ether_oxygens").html(physicochemical_properties["Ether Oxygens"]);
    $("#hydroxy_groups").html(physicochemical_properties["Hydroxy Groups"]);
    $("#carboxylic_acids").html(physicochemical_properties["Carboxylic Acids"]);

    $("#main_molecular_weight").html(physicochemical_properties["Molecular Weight"].toFixed(3));
    $("#main_hbd").html(physicochemical_properties["Hydrogen Bond Donors"]);
    $("#main_hba").html(physicochemical_properties["Hydrogen Bond Acceptors"]);
    $("#main_formal_charge").html(physicochemical_properties["Formal Charge"]);

}

function fill_external_sources(sources) {
    var url_hash = {
        "HMDB Accession" : "http://www.hmdb.ca/metabolites/",
        "PubChem ID" : "https://pubchem.ncbi.nlm.nih.gov/compound/",
        "CAS" : "http://www.molbase.com/en/cas-",
        "KEGG Compound" : "http://www.genome.jp/dbget-bin/www_bget?",
        "Wikidata" : "https://www.wikidata.org/wiki/",
        "Chemspider" : "http://www.chemspider.com/Chemical-Structure.",
        "BioCyc" : "http://www.biocyc.org/META/new-image?object="
    };

    for (var source in sources) {
        if (sources[source] != null) {
            var source_html = "<li class='list-group-item'>";
            source_html += source+": ";
            source_html += "<a href='";
            source_html += url_hash[source]+sources[source]+"' target='_blank'>";
            source_html += sources[source]+"</a></li>";
            $("#external_sources_list").append(source_html);
        }
    }
}

function supply_image(metabolite_id) {
    $("#structure").attr("src", getBaseURL()+"view/structure/"+metabolite_id);

}


function fill_isotopic_distribution_table(adduct_info) {
    $("#id_adduct_name").html(adduct_info["Type"]);
    $("#distribution_table tbody tr").remove();
    for (distribution_index in adduct_info["Isotopic Distribution"]) {
        var distribution = adduct_info["Isotopic Distribution"][distribution_index];
        $("#distribution_table").append("<tr><td class='text-center'>"+distribution[0].toFixed(3)+"</td><td class='text-center'>"+distribution[1].toFixed(2)+"</td></tr>")
    }

}

function chart_distribution(distribution_array) {
    var masses = [];
    var intensities = [];
    for (i in distribution_array) {
        masses.push(distribution_array[i][0]);
        intensities.push(distribution_array[i][1]);
    }


    var ctx = $("#myChart");

    var myChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels : masses,
            datasets : [{
                data : intensities
            }]
        },
        options: {
            responsive : true,
            maintainAspectRatio: false,
            legend : {
                display : false
            },
            scales: {
                xAxes: [{
                    display: true,
                    scaleLabel: {
                        display: true,
                        labelString: "Mass-to-ion (m/z)"
                    }
                }],
                yAxes: [{
                    display: true,
                    scaleLabel: {
                        display: true,
                        labelString: "Relative Intensity (%)"
                    }
                }]
            }
        }
    });

}

function fill_adducts_information(ionisation, adduct_information) {
    $("#ionisation_title").html(ionisation);

    $("#adduct_selector").empty();

    for (var adduct_index in adduct_information) {
        var adduct = adduct_information[adduct_index];
        var adduct_html = "<option ";
        adduct_html += "value='"+adduct_index+"'>";
        adduct_html += adduct["Type"] +": ";
        adduct_html += adduct["Accurate Mass"].toFixed(3) + " m/z";
        adduct_html += "</option>";
        $("#adduct_selector").append(adduct_html);
    }


}


function generate_chemical_formula_search_table(api_url) {
    $('#molecular_search_table').DataTable({
            "destroy": true,
            "searching": false,
            "lengthChange": false,
            "pageLength": 5,
            "ajax": {
                "url": encodeURI(api_url),
                "dataSrc": "_items"
            },
            "columns": [
                {
                    "title": "Metabolite Name",
                    "width": "40%",
                    "data": "Identification Information.Name",
                    "render": function (data, type, row) {
                        return "<a href='" + getBaseURL() + "view/" + row._id + "' target='_blank'>" + data + "</a>"
                    }
                },
                {
                    "title": "Molecular Formula",
                    "width": "10%",
                    "className": "dt-center",
                    "data": "Identification Information.Molecular Formula",
                    "render": function (data, type, row) {
                        return data.replace(/([0-9]+)/g, '<sub>$1</sub>');
                    }
                },
                {
                    "title": "Molecular Weight (g/mol)",
                    "data": "Physicochemical Properties",
                    "className": "dt-center",
                    "width": "10%",
                    "render": function (data, type, row) {
                        return data["Molecular Weight"].toFixed(3);
                    }
                }
            ],
        });
}


function fill_pathway_data(pathways) {



    function append_generator(pid, identifier) {
        var pathway_urls = {
            "KEGG" : "http://www.genome.jp/dbget-bin/www_bget?",
            "SMPDB" : "http://www.smpdb.ca/view/",
            "BioCyc" : "https://biocyc.org/META/NEW-IMAGE?type=PATHWAY&object="
        }

        var list_item = "<a href='"+pathway_urls[identifier]+pid+"'>";
        list_item += "<li class='list-group-item'>";
        if (identifier == "KEGG") {
            list_item += kegg_pathway_naming[pid];
        }

        else if (identifier == "SMPDB") {
            list_item += smpdb_pathway_naming[pid];
        }
        list_item += "</li></a>";
        return list_item
    }
    var kegg_pathways = pathways["KEGG"];
    var smpdb_pathways = pathways["SMPDB"];
    var biocyc_pathways = pathways["BioCyc"];

    if (kegg_pathways.length > 0) {
        $("#no_pathways").fadeOut(0);
        $(".kegg_pcount").text(kegg_pathways.length);
        $("#kegg_pathways").fadeIn(0);

        for (k_idx in kegg_pathways) {
            $("#kegg_list_group").append(append_generator(kegg_pathways[k_idx], "KEGG"))
        }
    }

    if (smpdb_pathways.length > 0) {
        $("#no_pathways").fadeOut(0);
        $(".smpdb_pcount").text(smpdb_pathways.length);
        $("#smpdb_pathways").fadeIn(0);

        for (smpdb_idx in smpdb_pathways) {
            $("#smpdb_list_group").append(append_generator(smpdb_pathways[smpdb_idx], "SMPDB"))
        }
    }

    if (biocyc_pathways.length > 0) {
        $("#no_pathways").fadeOut(0);
        $(".biocyc_pcount").text(biocyc_pathways.length);
        $("#biocyc_pathways").fadeIn(0);

        for (biocyc_idx in biocyc_pathways) {
            $("#biocyc_list_group").append(append_generator(biocyc_pathways[biocyc_idx], "BioCyc"))
        }
    }
}

function fill_tproperties_information(properties) {
    $("#hmdb_origins").html(array_to_text(properties["HMDB"]["Origins"]));
    $("#hmdb_bioloc").html(array_to_text(properties["HMDB"]["Biofluid Locations"]));
    $("#hmdb_tisloc").html(array_to_text(properties["HMDB"]["Tissue Locations"]));
}


function skeletons(inchikey) {
    var api_url = getBaseURL() + "api/metabolites/?where={";
    api_url += '"_id" : { "$regex" : "^' + inchikey.split("-")[0] +'"}}';
    api_url += '&projection={"Identification Information.Name" : 1}';

    var metabolites = get_metabolites(api_url);

    for (index in metabolites) {
        var m = metabolites[index];
        if (m["_id"] != inchikey) {
           $("#skeleton_list").append("<a href='" + getBaseURL() +"view/" + m["_id"] + "' target='_blank'><li class='list-group-item'>"
            +m["Identification Information"]["Name"] + "</li></a>")
        }


    }
}


function render_metabolite_view(metabolite_id) {
    var metabolite = get_metabolite(metabolite_id);
    fill_identification_infomation(metabolite["Identification Information"], metabolite_id);
    fill_physicochemical_properties(metabolite["Physicochemical Properties"]);
    fill_external_sources(metabolite["External Sources"]);
    supply_image(metabolite_id);

    fill_adducts_information("Neutral", metabolite["Adducts"]["Neutral"]);
    fill_isotopic_distribution_table(metabolite["Adducts"]["Neutral"][0]);
    chart_distribution(metabolite["Adducts"]["Neutral"][0]["Isotopic Distribution"]);
    fill_pathway_data(metabolite["Pathways"]);
    fill_tproperties_information(metabolite["Taxonomic Properties"]);

    skeletons(metabolite["_id"]);

    $("input[type=radio][name=ionisation]").change(function () {
        fill_adducts_information(this.value, metabolite["Adducts"][this.value]);
    });

    $("#adduct_selector").on("change", function () {
        var adduct_index = this.value;
        var ionisation = $("input[name=ionisation]:checked").val();
        fill_isotopic_distribution_table(metabolite["Adducts"][ionisation][adduct_index]);
        chart_distribution(metabolite["Adducts"][ionisation][adduct_index]["Isotopic Distribution"]);
    });

    $("#formula_search_button").click(function () {
        var api_url = getBaseURL() + "api/metabolites/?where={";
        api_url += '"Identification Information.Molecular Formula" : "';
        api_url += metabolite["Identification Information"]["Molecular Formula"] + '"}';
        api_url += '&projection={"Identification Information.Name" : 1, ';
        api_url += '"Identification Information.Molecular Formula" : 1, "Physicochemical Properties.Molecular Weight" : 1}&max_results=1000}';
        generate_chemical_formula_search_table(api_url);
        $('#formula_search_modal').modal('toggle');
    });
}

$("#spooky").click(function () {
    var spooky_array = [
        "thank mr skelator",
        "updoot 4 healthy bones",
        "spooky scary skeletons",
        "send shivers down your spine!"
    ];
    $("#spooky_quote").html(spooky_array[Math.floor(Math.random() * spooky_array.length)]);
    $("#spooky_skeleton").fadeIn(500).delay(1000).fadeOut(500);
});


$("#structure").on('click', function() {
    $('.enlargeImageModalSource').attr('src', $(this).attr('src'));
	$('#enlargeImageModal').modal('show');
});
