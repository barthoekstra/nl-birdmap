-- Amsterdam migration project 
hostvar="robin1.e-ecology.nl"
uservar="johannes_groeve"
pwvar=""

PGPASSWORD=$pwvar psql -h $hostvar -U $uservar -d vp_db -c "
DROP TABLE temp.ams_migration_altitude_profiles;
CREATE TABLE temp.ams_migration_altitude_profiles AS
WITH x AS 
(
	SELECT *, unnest(vp_file_id_array) vp_file_id 
	FROM temp.ams_migration_hmd 
)
SELECT vp_versioning_id, date_trunc, odim_code, year_season, year, season, day_night, mtr, prop, mtr_cum, prop_cum, top, n_vp,
hght,
sum(dens::numeric) dens,
sum(n) n,
avg(dens::numeric) dens_avg,
avg(n) n_avg
FROM x 
JOIN main.vp_data USING (vp_file_id) 
GROUP BY vp_versioning_id, date_trunc, odim_code, year_season, year, season, day_night, hght, mtr, prop, mtr_cum, prop_cum, top, n_vp
ORDER BY vp_versioning_id, odim_code, year_season, year, season, day_night, top 
;"


COMMENT ON TABLE temp.ams_migration_altitude_profiles IS $$

For every full day and night period the total and average density (dens, dens_avg) and number of birds (n, n_avg) per altitude bin was calculated. 
These were calculated by joining the MTR ordered seasonal day/night migration dataset generated in [step 1](https://ams-migration-uva-ibed-ame-multisensor-050112b4bde163c8deefb305.gitlab.io/vp_mtr_ordered_seasonal_day_night_migration.html) 
with the vp_data using the vertical profile identifiers (vp_file_id) stored as an array in the schema.table temp.ams_migration_hmd and by aggregating using the columns odim_code, date_trunc, season, year and day_night. The daily altitude profiles can be 
used to filter for the relevant high migration days (cumulative MTR < 0.50) and further aggregate at different time scales (e.g. season, year). 
$$;

COMMENT ON COLUMN temp.ams_migration_altitude_profiles.vp_versioning_id IS 'which versioning is used to compute HMD: default (0) or h50layer (28)';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.date_trunc IS 'reference night / day. The evening is assigned to the date of the next night';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.odim_code IS 'radar odim code';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.year_season IS 'year and season';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.year IS 'year';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.season IS 'season';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.day_night IS 'whether the MTR is calculated for the Day or Night period as defined from sunrise / sunset';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.mtr IS 'total MTR for the day or night period (between sunrise and sunset and vice versa)';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.prop IS 'proportion of the total MTR for the day and night contributing to the full year-season';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.mtr_cum IS 'cumulative total MTR for the day and night period';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.prop_cum IS 'cumulative proportion of the total MTR for the day and night period contributing to the full year-season';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.top IS 'total MTR based index ordered by descreasing migration intensity per year-season';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.n_vp IS 'number of vertical profiles within the time range';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.dens IS 'MTR sum for the day and night period per altitude bin';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.n IS 'counts sum for the day and night period per altitude bin';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.dens_avg IS 'MTR average for the day and night period per altitude bin';
COMMENT ON COLUMN temp.ams_migration_altitude_profiles.n_avg IS 'counts average for the day and night period per altitude bin';



	