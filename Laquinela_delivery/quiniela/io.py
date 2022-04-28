import sqlite3

import pandas as pd

import settings


def load_matchday(season, division, matchday):
    season = season.split('-')[0]
    with sqlite3.connect(settings.DATABASE_PATH) as conn:
        data = pd.read_sql(f"""
            SELECT * FROM Matches
                WHERE season = '{season}'
                  AND division = {division}
                  AND matchday = {matchday}
        """, conn)
    if data.empty:
        raise ValueError("There is no matchday data for the values given")
    return data


def load_historical_data(seasons):
    with sqlite3.connect(settings.DATABASE_PATH) as conn:
        if seasons == "all":
            data = pd.read_sql("SELECT * FROM Matches", conn)
        else:
            seasons = [season.split('-')[0] for season in seasons]
            data = pd.read_sql(f"""
                SELECT * FROM Matches
                    WHERE season IN {tuple(seasons)}
            """, conn)
    if data.empty:
        raise ValueError(f"No data for seasons {seasons}")
    return data


def save_predictions(predictions):
    cols_to_save = ['season', 'division', 'matchday', 'home_team', 'score', 'away_team', 'pred']
    with sqlite3.connect(settings.DATABASE_PATH) as conn:
        predictions[cols_to_save].to_sql(name="Predictions", con=conn, if_exists="append", index=False)
