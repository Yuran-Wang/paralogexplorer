import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from shiny import App, Inputs, Outputs, Session, reactive, ui, render
from shinywidgets import output_widget, render_widget
from taigapy import create_taiga_client_v3

# Initialize Taiga client
taiga_client = create_taiga_client_v3()

# Load Taiga Dataset
paralog_guide_map = taiga_client.get(name='paralogscreen06202024-816a', version=29, file='ParalogFullGuideMap')
paralog_screen_map = taiga_client.get(name='paralogscreen06202024-816a', version=29, file='ParalogScreenMap')
paralog_gene_effect = taiga_client.get(name='paralogscreen06202024-816a', version=29, file='ParalogGeneEffect')
paralog_lfc_gene_screen = taiga_client.get(name='paralogscreen06202024-816a', version=29, file='ParalogFullLfcGeneScreen')
paralog_lfc_gene_screen.set_index('GuideTargetSymbol', inplace=True)
paralog_lfc_sgrna_screen = taiga_client.get(name='paralogscreen06272024-fa85', version=9, file='ParalogFullLfcSgrnaScreen')

PRMT5_group = '/Users/wangyura/Downloads/combined_sheets/PRMT5_combined.csv'
PT2385_group = '/Users/wangyura/Downloads/PT2385_combined.csv'
ICC137_group = '/Users/wangyura/Downloads/combined_sheets/ICC137_combined.csv'
parp_parg_group = '/Users/wangyura/Downloads/combined_sheets/parp_parg_combined.csv'
ICC21_group = '/Users/wangyura/Downloads/combined_sheets/ICC21_combined.csv'
BI_group = '/Users/wangyura/Downloads/combined_sheets/BI_combined.csv'

# Extract 'Anchor' type rows
anchor_data = paralog_screen_map[paralog_screen_map['ScreenType'] == 'Anchor']

# Filter ParalogFullGuideMap to include only rows with "_" in GuideTargetSymbol
filtered_guide_map = paralog_guide_map[paralog_guide_map['GuideTargetSymbol'].str.contains('_')]

# Define UI
app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.card(
                ui.h2("Paralog Screen Selector"),
                ui.input_selectize(
                    "group_type",
                    "Select group",
                    choices={"PRMT5": "PRMT5", "ICC137": "ICC137", "PT2385": "PT2385", "parp_parg": "parp_parg", "ICC21": "ICC21", "BI": "BI"},
                    multiple=True,
                    selected=["PRMT5"]
                ),
                ui.output_data_frame("comparison_table"),
                ui.output_text_verbatim("all_comparisons"),     
            ),
            ui.card(
                ui.card_header("Box Plot Gene Selection"),
                ui.output_data_frame("gene_selection_table"),
                ui.output_text_verbatim("selected_rows_box_gene"),
                ui.output_text_verbatim("display_selected_rows"),  # New output to display selected rows
                class_="mt-4"
            )
        ),
        ui.panel_main(
            ui.accordion(
                ui.accordion_panel(
                    "Waterfall Plot",
                    ui.card(
                        ui.card_header(
                            "Waterfall Plot",
                            ui.popover(
                                "Options",
                                ui.input_radio_buttons(
                                    "waterfall_plot_data",
                                    "Choose calculation:",
                                    ["ParalogFullLfc", "Beta"],
                                    selected="Beta",
                                    inline=True,
                                ),
                                title="Choose calculation",
                            ),
                            class_="d-flex justify-content-between align-items-center",
                        ),
                        output_widget("waterfall_plot"),
                        full_screen=True,
                    ),
                ),
                ui.accordion_panel(
                    "Box Plot",
                    ui.card(
                        ui.card_header(
                            "Box Plot",
                            ui.popover(
                                "Options",
                                ui.input_radio_buttons(
                                    "box_plot_gene",
                                    "Select Box Plot number:",
                                    ["1", "2"],
                                    selected="1",
                                    inline=True,
                                ),
                                title="Choose gene",
                            ),
                            class_="d-flex justify-content-between align-items-center",
                        ),
                        output_widget("box_plot"),
                        full_screen=True,
                    ),
                ),
                title="Paralog Visualization"
            )
        )
    )
)

# Define server logic
def server(input: Inputs, output: Outputs, session: Session):

    # Define reactive calculations
    @reactive.Calc
    def selected_rows():
        return gene_selection_table.cell_selection()["rows"]

    @reactive.Calc
    def selected_data():
        rows = selected_rows()
        if not rows:
            return pd.DataFrame()
        return filtered_guide_map.iloc[list(rows)]

    @reactive.Calc
    def guide_labels_combined():
        selected = selected_data()
        guide_labels = []
        for symbol in selected['GuideTargetSymbol'].unique():
            # Search for the combined symbol
            guide_labels += paralog_guide_map[paralog_guide_map['GuideTargetSymbol'] == symbol]['GuideLabel'].tolist()
        return guide_labels

    @reactive.Calc
    def guide_labels_separate1():
        selected = selected_data()
        guide_labels = []
        for symbol in selected['GuideTargetSymbol'].unique():
            # Split the symbol and search for the first part before "_"
            part1 = symbol.split('_')[0]
            guide_labels += paralog_guide_map[paralog_guide_map['GuideTargetSymbol'] == part1]['GuideLabel'].tolist()
        return guide_labels

    @reactive.Calc
    def guide_labels_separate2():
        selected = selected_data()
        guide_labels = []
        for symbol in selected['GuideTargetSymbol'].unique():
            # Split the symbol and search for the second part after "_"
            parts = symbol.split('_')
            if len(parts) > 1:
                part2 = parts[1]
                guide_labels += paralog_guide_map[paralog_guide_map['GuideTargetSymbol'] == part2]['GuideLabel'].tolist()
        return guide_labels

    @reactive.Calc
    def full_identifier_values():
        combined_labels = guide_labels_combined()
        separate1_labels = guide_labels_separate1()
        separate2_labels = guide_labels_separate2()
        all_labels = combined_labels + separate1_labels + separate2_labels
        return [value.strip().replace('"', '') for value in ', '.join(all_labels).split(',') if value.strip()]

    @reactive.Calc
    def transformed_gene_df():
        return pd.DataFrame(full_identifier_values(), columns=["GuideLabel"])
    
    # Function to check if cell_line_column exists
    def check_column_existence(df, column_name):
        return column_name in df.columns

    @reactive.Calc
    def transformed_df():
        group_types = input.group_type()
        dfs = []
        for group_type in group_types:
            if group_type == "PRMT5":
                df = pd.read_csv(PRMT5_group)
            elif group_type == "ICC137":
                df = pd.read_csv(ICC137_group)
            elif group_type == "PT2385":
                df = pd.read_csv(PT2385_group)
            elif group_type == "parp_parg":
                df = pd.read_csv(parp_parg_group)
            elif group_type == "ICC21":
                df = pd.read_csv(ICC21_group)
            elif group_type == "BI":
                df = pd.read_csv(BI_group)
            else:
                df = pd.DataFrame()
            dfs.append(df)

        df = pd.concat(dfs, ignore_index=True)
        unique_experiments = df["Experiment"].unique()
        selected_data = pd.DataFrame(unique_experiments, columns=["Experiment"])

        lines = selected_data.to_string(index=False).split('\n')
        filtered_lines = [line.strip() for line in lines if "Experiment" not in line.strip()]

        transformed_text = []
        for line in filtered_lines:
            parts = line.split('_vs_')
            if len(parts) == 2:
                base, drug_suffix = parts[0].rsplit('_', 1)
                dmso_suffix = parts[1]
                transformed_text.append([f"{base}.{drug_suffix}", base, drug_suffix])
                transformed_text.append([f"{base}.{dmso_suffix}", base, dmso_suffix])

        return pd.DataFrame(transformed_text, columns=["Screen", "CellLine", "Treatment"])
    
    @reactive.Calc
    def selected_calculation():
        return input.waterfall_plot_data()
    
    # Set the frame for filterable and selectable data, ensuring unique values
    @output
    @render.data_frame
    def gene_selection_table():
        unique_guide_map = filtered_guide_map[['GuideTargetSymbol']].drop_duplicates()
        return render.DataGrid(unique_guide_map, filters=True, selection_mode="rows")


    # Comparison Table without row selection
    @output
    @render.data_frame
    def comparison_table():
        group_types = input.group_type()
        dfs = []
        for group_type in group_types:
            if group_type == "PRMT5":
                df = pd.read_csv(PRMT5_group)
            elif group_type == "ICC137":
                df = pd.read_csv(ICC137_group)
            elif group_type == "PT2385":
                df = pd.read_csv(PT2385_group)
            elif group_type == "parp_parg":
                df = pd.read_csv(parp_parg_group)
            elif group_type == "ICC21":
                df = pd.read_csv(ICC21_group)
            elif group_type == "BI":
                df = pd.read_csv(BI_group)
            else:
                df = pd.DataFrame()
            dfs.append(df)

        df = pd.concat(dfs, ignore_index=True)
        unique_comparisons = df["Experiment"].unique()
        return render.DataGrid(pd.DataFrame(unique_comparisons, columns=["Experiment"]), filters=True, height='400px')

    @output
    @render.text
    def all_comparisons():
        transformed_data = transformed_df()
        print("text verbatim for comparison:", transformed_data)
        return transformed_data.to_string(index=False)

    @output
    @render.text
    def selected_rows_box_gene():
        rows = selected_rows()
        if not rows:
            return "No rows selected."
        else:
            data = selected_data()

            # Initialize dictionaries to store the labels
            combined_guide_labels = {}
            separate1_guide_labels = {}
            separate2_guide_labels = {}

            for symbol in data['GuideTargetSymbol'].unique():
                combined_guide_labels[symbol] = paralog_guide_map[paralog_guide_map['GuideTargetSymbol'] == symbol]['GuideLabel'].tolist()
                
                # Split the symbol into parts
                parts = symbol.split('_')
                first_part = parts[0]
                second_part = parts[1] if len(parts) > 1 else ''
                
                # Print the full symbol, first part, and second part
                print("Full GuideTargetSymbol:", symbol)
                print("First part:", first_part)
                print("Second part:", second_part)
                
                # Get separate guide labels (first and second parts)
                separate1_guide_labels[first_part] = paralog_guide_map[paralog_guide_map['GuideTargetSymbol'] == first_part]['GuideLabel'].tolist()
                if second_part:
                    separate2_guide_labels[second_part] = paralog_guide_map[paralog_guide_map['GuideTargetSymbol'] == second_part]['GuideLabel'].tolist()

            # Function to format lists as strings
            def format_guide_labels(label_dict):
                formatted_str = ""
                for key, values in label_dict.items():
                    cleaned_values = [value.strip().replace('"', '') for value in ', '.join(values).split(',') if value.strip()]
                    formatted_list = '[\n    ' + ',\n    '.join(f'"{value}"' for value in cleaned_values) + '\n]'
                    formatted_str += f"{key}:\n{formatted_list}\n\n"
                return formatted_str

            # Format the output
            combined_formatted_str = format_guide_labels(combined_guide_labels)
            separate1_formatted_str = format_guide_labels(separate1_guide_labels)
            separate2_formatted_str = format_guide_labels(separate2_guide_labels)

            # Return combined and separate lists as formatted strings
            return f"Combined Guide Labels:\n{combined_formatted_str}\nSeparate Guide Labels 1:\n{separate1_formatted_str}\nSeparate Guide Labels 2:\n{separate2_formatted_str}"

    
    @output
    @render_widget
    def waterfall_plot():
        outputdata = transformed_df()
        print("Waterfall data:", outputdata)
        if outputdata.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No data found for the selection",
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=20)
            )
            return fig

        group_types = input.group_type()
        dfs = []

        for group_type in group_types:
            if group_type == "PRMT5":
                df = pd.read_csv(PRMT5_group)
            elif group_type == "ICC137":
                df = pd.read_csv(ICC137_group)
            elif group_type == "PT2385":
                df = pd.read_csv(PT2385_group)
            elif group_type == "parp_parg":
                df = pd.read_csv(parp_parg_group)
            elif group_type == "ICC21":
                df = pd.read_csv(ICC21_group)
            elif group_type == "BI":
                df = pd.read_csv(BI_group)
            else:
                df = pd.DataFrame()
            dfs.append(df)

        # Filter and prepare the list of dataframes for plotting
        df_list = []

        for df in dfs:
            if df.empty:
                continue

            # Get unique combinations of Treatment and cellLine
            unique_combinations = df[['Treatment', 'cellLine']].drop_duplicates()

            for _, row in unique_combinations.iterrows():
                treatment = row['Treatment']
                cell_line = row['cellLine']
                filtered_data = df[(df['Treatment'] == treatment) & (df['cellLine'] == cell_line)]
                
                if not filtered_data.empty:
                    df_list.append(filtered_data)

        # Function to create the water flow plot
        def water_flow_plot_v2(df_list, calculation, title=None, width=800, height=800):
            if len(df_list) == 0:
                fig = go.Figure()
                fig.add_annotation(
                    text="No data available for waterfall plot",
                    xref="paper", yref="paper",
                    showarrow=False,
                    font=dict(size=20)
                )
                return fig

            # Generate the waterfall plot
            fig = make_subplots(
                rows=len(df_list),
                cols=2,
                subplot_titles=[f'{df["Treatment"].iloc[0]} - {df["cellLine"].iloc[0]}' for df in df_list for _ in range(2)]
            )

            for idx, df in enumerate(df_list):
                print(f"Creating plot for: {df['Treatment'].iloc[0]} - {df['cellLine'].iloc[0]}")
                print(f"Number of rows: {len(df)}")
                
                # Check if the required columns are present
                if calculation == "Beta":
                    x_col = 'rank_beta'
                    y_col = 'Treatment|beta'
                else:  # ParalogFullLfc
                    x_col = 'rank_LFC_diff'
                    y_col = 'LFC_diff'

                hover_indices = [0, 1, 2, 3, 4, 6, 7]
                hover_columns = [df.columns[i] for i in hover_indices if i < len(df.columns)]
                
                # Paralog Plot
                paralog_df = df[df['gene_type'] == 'paralog']
                print(f"Number of Paralog rows: {len(paralog_df)}")
                if not paralog_df.empty:
                    for trace in px.scatter(
                        paralog_df,
                        x=x_col,
                        y=y_col,
                        color_discrete_sequence=['blue'],
                        hover_data=hover_columns
                    ).data:
                        fig.add_trace(trace, row=idx + 1, col=1)

                # SingleGene Plot
                singlegene_df = df[df['gene_type'] == 'single']
                print(f"Number of SingleGene rows: {len(singlegene_df)}")
                if not singlegene_df.empty:
                    for trace in px.scatter(
                        singlegene_df,
                        x=x_col,
                        y=y_col,
                        color_discrete_sequence=['red'],
                        hover_data=hover_columns
                    ).data:
                        fig.add_trace(trace, row=idx + 1, col=2)

            fig.update_layout(height=400 * len(df_list), showlegend=True, title_text=title)
            return fig

        selected_calc = selected_calculation()

        # Call the function with the list of dataframes and selected calculation
        fig = water_flow_plot_v2(df_list, calculation=selected_calc)
        return fig



    @output
    @render_widget
    def box_plot():
        screen_columns = transformed_df()['Screen'].tolist()
        unique_cell_lines = transformed_df()['CellLine'].unique()
        num_unique_cell_lines = len(unique_cell_lines)
        print("Unique Cell Lines: ", unique_cell_lines)


         # Verify the existing columns
        print("Original columns:", paralog_lfc_sgrna_screen.columns)

        # Define a dictionary with the old and new column names
        rename_columns = {
            'BI_SW1573.AMG510': 'SW1573.AMG510',
            'BI_SW1573.DMSO': 'SW1573.DMSO.BI',
            'RCC10RGB.PT2385': 'RCC10RGB.PT',
            'BI_SW1573.BI3406': 'SW1573.BI3406',
            'BI_SW1573.Combo': 'SW1573.Combo'
        }

        # Create a list of columns to use, renaming if necessary
        comparison_columns = []
        for cell_line in unique_cell_lines:
            if check_column_existence(paralog_lfc_sgrna_screen, cell_line):
                comparison_columns.append(cell_line)

        # Rename the screen columns based on the dictionary if they exist
        renamed_screen_columns = []
        for col in screen_columns:
            renamed_col = rename_columns.get(col, col)
            if renamed_col in paralog_lfc_sgrna_screen.columns:
                renamed_screen_columns.append(renamed_col)
            elif col in paralog_lfc_sgrna_screen.columns:
                renamed_screen_columns.append(col)

        comparison_columns += renamed_screen_columns

        print("Columns to use:", comparison_columns)
        print("Screen used:", screen_columns)
        print("Unique cell lines:", unique_cell_lines)
        print('Treatment used:', transformed_df()['Treatment'].unique())


        combined_identifier_values = guide_labels_combined()
        separate1_identifier_values = guide_labels_separate1()
        separate2_identifier_values = guide_labels_separate2()

        combined_matching_indices = [value for value in combined_identifier_values if value in paralog_lfc_sgrna_screen.index]
        print("Combined Matching indices: ", combined_matching_indices)
        separate1_matching_indices = [value for value in separate1_identifier_values if value in paralog_lfc_sgrna_screen.index]
        print("Separate 1 indicies: ", separate1_matching_indices)
        separate2_matching_indices = [value for value in separate2_identifier_values if value in paralog_lfc_sgrna_screen.index]
        print("Separate 2 indicies: ", separate2_matching_indices)

        if not combined_matching_indices and not separate1_matching_indices and not separate2_matching_indices:
            fig = go.Figure()
            fig.add_annotation(
                text="No data found for the gene",
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=20)
            )
            return fig

        # Define the symbols based on the selected data
        selected = selected_data()
        if 'GuideTargetSymbol' in selected.columns:
            symbol = selected['GuideTargetSymbol'].iloc[0]
        else:
            # Handle the case where 'GuideTargetSymbol' column is missing
            symbol = "Unknown_Unknown"
            
        first_part = symbol.split('_')[0]
        second_part = symbol.split('_')[1] if '_' in symbol else ''

        guide_label_dict = {value: symbol for value in combined_identifier_values}
        guide_label_dict.update({value: first_part for value in separate1_identifier_values})
        guide_label_dict.update({value: second_part for value in separate2_identifier_values})

        combined_matching_rows = paralog_lfc_sgrna_screen.loc[combined_matching_indices, comparison_columns] if combined_matching_indices else pd.DataFrame()
        separate1_matching_rows = paralog_lfc_sgrna_screen.loc[separate1_matching_indices, comparison_columns] if separate1_matching_indices else pd.DataFrame()
        separate2_matching_rows = paralog_lfc_sgrna_screen.loc[separate2_matching_indices, comparison_columns] if separate2_matching_indices else pd.DataFrame()

        combined_matching_rows = combined_matching_rows.reset_index().melt(id_vars="index", value_vars=comparison_columns, var_name="Condition", value_name="Value")
        separate1_matching_rows = separate1_matching_rows.reset_index().melt(id_vars="index", value_vars=comparison_columns, var_name="Condition", value_name="Value")
        separate2_matching_rows = separate2_matching_rows.reset_index().melt(id_vars="index", value_vars=comparison_columns, var_name="Condition", value_name="Value")

        combined_matching_rows['Type'] = combined_matching_rows['Condition'].apply(lambda x: next((treatment for treatment in transformed_df()['Treatment'].unique() if treatment in x), 'Unknown'))
        separate1_matching_rows['Type'] = separate1_matching_rows['Condition'].apply(lambda x: next((treatment for treatment in transformed_df()['Treatment'].unique() if treatment in x), 'Unknown'))
        separate2_matching_rows['Type'] = separate2_matching_rows['Condition'].apply(lambda x: next((treatment for treatment in transformed_df()['Treatment'].unique() if treatment in x), 'Unknown'))

        combined_matching_rows['GuideLabel'] = combined_matching_rows['index'].map(guide_label_dict)
        separate1_matching_rows['GuideLabel'] = separate1_matching_rows['index'].map(guide_label_dict)
        separate2_matching_rows['GuideLabel'] = separate2_matching_rows['index'].map(guide_label_dict)

        combined_matching_rows['Condition'] = symbol + '-' + combined_matching_rows['Condition']
        separate1_matching_rows['Condition'] = first_part + '-' + separate1_matching_rows['Condition']
        separate2_matching_rows['Condition'] = second_part + '-' + separate2_matching_rows['Condition']

        all_matching_rows = pd.concat([separate1_matching_rows, combined_matching_rows, separate2_matching_rows])
      # Add the "CellLine" column with a special case for values containing 'SW1573'
        all_matching_rows['CellLine'] = all_matching_rows['Condition'].apply(
            lambda x: 'BI_SW1573' if 'SW1573.DMSO.BI' in x else
                    ('BI_SW1573' if 'SW1573.BI3406' in x else
                    ('BI_SW1573' if 'SW1573.Combo' in x else 
                    ('BI_SW1573' if 'SW1573.AMG510' in x else 
                    (x.split('-')[1].split('.')[0] if '-' in x else 'N/A'))))
        )

        all_matching_rows.sort_values(by=['CellLine', 'Type', 'GuideLabel'], inplace=True)
        unique_guidelabels = all_matching_rows['GuideLabel'].unique()

        # Add the "Guide Target Symbol" column
        all_matching_rows['Guide Target Symbol'] = all_matching_rows['index'].map(guide_label_dict)
        # Set display options to show the entire DataFrame
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_colwidth', None)
        pd.set_option('display.width', None)
        print("dataframe: All Matching rows: ", all_matching_rows)

        # Map treatment types to colors
        color_map = {
            'DMSO': 'lightgrey',
            'MRTX9768': 'red',
            'MRTX1719': 'purple',
            'AG43192': 'orange',
            'Unknown': 'green',
            '1nM.RLY4008': 'pink',
            '2.5nM.RLY4008': 'yellow',
            'Niraparib': 'peru',
            'PARGi': 'skyblue',
            'PARP': 'violet',
            'Pemigatinib': 'lightsalmon',
            'PARG': 'limegreen',
            'AMG510': 'deepskyblue',
            'BI3406': 'orchid',
            'COMBO': 'gold',
            'Combo': 'turquoise',
            'PT2385': 'paleturquoise'

        }

        # Create subplots with all facets horizontally next to each other
        fig = make_subplots(rows=1, cols=num_unique_cell_lines, subplot_titles=unique_cell_lines)

        min_y = float('inf')
        max_y = float('-inf')
        box_width = 0.5

        # Initialize a set to keep track of added legend items
        added_legend_items = set()

        for i, cell_line in enumerate(unique_cell_lines):
            df_plot = all_matching_rows[all_matching_rows['CellLine'] == cell_line]
            for guide_label in df_plot['GuideLabel'].unique():
                sub_df = df_plot[df_plot['GuideLabel'] == guide_label]
                for treatment in sub_df['Type'].unique():
                    treatment_df = sub_df[sub_df['Type'] == treatment]
                    min_y = min(min_y, treatment_df['Value'].min())
                    max_y = max(max_y, treatment_df['Value'].max())

        for i, cell_line in enumerate(unique_cell_lines):
            df_plot = all_matching_rows[all_matching_rows['CellLine'] == cell_line]
            x_ticks = []
            tick_labels = []
            guide_label_centers = {}

            # Create distinct x-axis categories for each subplot
            for guide_label in df_plot['GuideLabel'].unique():
                sub_df = df_plot[df_plot['GuideLabel'] == guide_label]
                treatment_positions = []
                for treatment in sub_df['Type'].unique():
                    treatment_df = sub_df[sub_df['Type'] == treatment]
                    # Create a unique tick value for each combination
                    tick_val = f"{guide_label}_{treatment}"
                    x_ticks.append(tick_val)
                    treatment_positions.append(tick_val)
                    fig.add_trace(
                        go.Box(
                            x=[tick_val]*len(treatment_df),
                            y=treatment_df['Value'],
                            name=f"{guide_label} - {treatment}",
                            marker=dict(color=color_map[treatment]),
                            boxpoints='all',  # Add dots
                            jitter=0.5, # Increase jitter to spread out the points
                            pointpos=0,
                            fillcolor='rgba(255,255,255,0)',
                            width=box_width, # Increase the width of the boxes
                            showlegend=False if treatment in added_legend_items else True,
                            legendgroup=treatment
                        ),
                        row=1,
                        col=i+1
                    )
                    if treatment not in added_legend_items:
                        added_legend_items.add(treatment)
                        # Add dummy trace for legend item
                        fig.add_trace(
                            go.Scatter(
                                x=[None],
                                y=[None],
                                mode='markers',
                                marker=dict(color=color_map[treatment]),
                                legendgroup=treatment,
                                showlegend=True,
                                name=treatment
                            )
                        )

                    # Update min_y and max_y based on data
                    min_y = min(min_y, treatment_df['Value'].min())
                    max_y = max(max_y, treatment_df['Value'].max())
                # Place the guide label in the middle of its treatment positions
                if treatment_positions:
                    center_pos = treatment_positions[len(treatment_positions) // 2]
                    guide_label_centers[guide_label] = center_pos
                tick_labels.extend([""] * len(treatment_positions))

            # Update x-axes to show only guide labels at their centers
            fig.update_xaxes(
                tickvals=x_ticks,
                ticktext=tick_labels,
                tickangle=45,
                row=1,
                col=i+1
            )

            # Add guide labels at their calculated center positions
            for guide_label, center_pos in guide_label_centers.items():
                fig.add_annotation(
                    x=center_pos,
                    y=min_y - 0.2,
                    text=guide_label,
                    showarrow=False,
                    xanchor='center',
                    yanchor='top',
                    font=dict(size=10),
                    xref=f'x{i+1}'
                )

        # Adjust the y-axis range to zoom in on the box plots
        fig.update_yaxes(range=[min_y - 0.2, max_y + 0.2])
        fig.add_hline(y=0, line_dash="dash", line_color="gray")
        fig.add_hline(y=-1, line_dash="dash", line_color="red")
        fig.update_layout(
            boxmode="group",
            margin={"l": 0, "r": 0, "t": 20, "b": 40},
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(size=10),
            legend_title='Treatment'
        )

        return fig



    @output
    @render.ui
    def plot_output():
        plot_type = input.plot_type()
        if plot_type == "Waterfall":
            return output_widget("waterfall_plot")
        elif plot_type == "Box Plot":
            return output_widget("box_plot")

app = App(app_ui, server)

if __name__ == "__main__":
    app.run()