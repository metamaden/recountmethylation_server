""" celeryconfig.py
    Configurations for celery task queue management. Specifies backend db.
"""

CELERY_RESULT_BACKEND = "database"
CELERY_RESULT_DBURI = "sqlite:///serverout.db"