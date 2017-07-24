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

function fill_physiochemical_properties(physiochemical_properties) {
    $("#molecular_weight").html(physiochemical_properties["Molecular Weight"].toFixed(3));
    $("#hbd").html(physiochemical_properties["Hydrogen Bond Donors"]);
    $("#hba").html(physiochemical_properties["Hydrogen Bond Acceptors"]);
    $("#rotatable_bonds").html(physiochemical_properties["Rotatable Bonds"]);
    $("#formal_charge").html(physiochemical_properties["Formal Charge"]);
    $("#fraction_of_sp3").html(physiochemical_properties["Fraction of SP3 Carbon"].toFixed(4));
    $("#clogp").html(physiochemical_properties["clogP"].toFixed(4));
    $("#polar_surface_area").html(physiochemical_properties["Polar Surface Area"].toFixed(4));
    $("#mr_values").html(physiochemical_properties["MR Values"].toFixed(4));
    $("#rings").html(physiochemical_properties["Rings"]);
    $("#aromatic_rings").html(physiochemical_properties["Aromatic Rings"]);
    $("#heavy_atoms").html(physiochemical_properties["Heavy Atoms"]);
    $("#secondary_amines").html(physiochemical_properties["Secondary Amines"]);
    $("#ether_oxygens").html(physiochemical_properties["Ether Oxygens"]);
    $("#hydroxy_groups").html(physiochemical_properties["Hydroxy Groups"]);
    $("#carboxylic_acids").html(physiochemical_properties["Carboxylic Acids"]);

    $("#main_molecular_weight").html(physiochemical_properties["Molecular Weight"].toFixed(3));
    $("#main_hbd").html(physiochemical_properties["Hydrogen Bond Donors"]);
    $("#main_hba").html(physiochemical_properties["Hydrogen Bond Acceptors"]);
    $("#main_formal_charge").html(physiochemical_properties["Formal Charge"]);

}

function fill_external_sources(sources) {
    var url_hash = {
        "HMDB Accession" : "http://www.hmdb.ca/metabolites/",
        "PubChem ID" : "https://pubchem.ncbi.nlm.nih.gov/compound/"
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
    $.ajax({
        url: getBaseURL() + "gen_structure/" + metabolite_id,
        success: function (result) {
            $("#structure").attr("src", "data:image/png;base64," + result);
        }
    });
}


function fill_isotopic_distribution_table(adduct_info) {
    $("#id_adduct_name").html(adduct_info["Type"]);
    $("#distribution_table tbody tr").remove();
    for (distribution_index in adduct_info["Isotopic Distribution"]) {
        var distribution = adduct_info["Isotopic Distribution"][distribution_index];
        $("#distribution_table").append("<tr><td class='text-center'>"+distribution[0].toFixed(3)+"</td><td class='text-center'>"+distribution[1].toFixed(2)+"</td></tr>")
    }

}

function fill_adducts_information(ionisation, adduct_information) {
    $("#ionisation_title").html(ionisation);
    $("#len_adducts").html(adduct_information.length);


    $("#adduct_list").empty();
    for (var adduct_index in adduct_information) {
        var adduct = adduct_information[adduct_index];
        var adduct_html = "<li class='list-group-item' id='adduct_select'";
        adduct_html += "name='"+adduct_index+"'";
        adduct_html +="><b>";
        adduct_html += adduct["Type"] +":</b> ";
        adduct_html += adduct["Accurate Mass"].toFixed(3);
        adduct_html += "</li>";
        $("#adduct_list").append(adduct_html);
    }


    console.log(adduct_information);

    $('#adduct_list').on('click', "#adduct_select", function() {
        var adduct_index = $(this).attr("name");
        fill_isotopic_distribution_table(adduct_information[adduct_index]);
    });
}


function render_metabolite_view(metabolite_id) {
    var metabolite = get_metabolite(metabolite_id);
    fill_identification_infomation(metabolite["Identification Information"], metabolite_id);
    fill_physiochemical_properties(metabolite["Physiochemical Properties"]);
    fill_external_sources(metabolite["External Sources"]);
    supply_image(metabolite_id);

    fill_adducts_information("Neutral", metabolite["Adducts"]["Neutral"]);
    fill_isotopic_distribution_table(metabolite["Adducts"]["Neutral"][0]);
    $("input[type=radio][name=ionisation]").change(function () {
        fill_adducts_information(this.value, metabolite["Adducts"][this.value])
    });



}


